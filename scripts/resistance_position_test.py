# posegate/scripts/resistance_position_test.py
"""Tests whether the miner's conserved-contact residues, collapsed onto
HIV-1 protease's mature-monomer numbering (1-99) and symmetry-collapsed
across the homodimer's two chains, are enriched for Stanford HIVdb's
known major protease-inhibitor resistance mutation positions.

This is a validation exercise, not a general-purpose feature: HIV
protease already has an authoritative, hand-curated resistance mutation
list (Stanford HIVdb), so it's the one target where "does a conserved
ligand-contact residue predict a resistance-relevant residue" can be
checked against real ground truth rather than asserted. If the
enrichment holds here, that's evidence the same automatic approach could
flag resistance-relevant residues on a target that does NOT have a
Stanford-HIVdb-equivalent curated resource.

Ensemble (14 structures, verified via RCSB): wild-type, ligand-bound HIV-1
protease only -- individually checked by RCSB title text to exclude
resistance mutants, apo structures, and other Gag-Pol domains (RNase H,
integrase) that share the same UniProt accession (P04585 is the whole
Gag-Pol polyprotein, not protease alone) but are a different protein
region entirely.

Numbering: SIFTS-remapping onto P04585 lands each residue at UniProt
position (mature-protease author number + 488) -- verified as a constant
offset across three independently-checked structures (1HVR, 1DMP, 5KAO;
see conversation). Converting back with -488 recovers Stanford's
1-99 mature-enzyme convention. Both chains of the homodimer map to the
SAME UniProt range (489-587), so subtracting 488 and dropping the chain
letter is exactly the symmetry-collapsing this analysis needs -- chain A
residue i and chain B residue i are the same physical monomer position.

Ground truth: MAJOR_RESISTANCE_POSITIONS is the Stanford HIVdb major PI
resistance mutation position list, confirmed via Shafer, "HIV-1 Protease
and Reverse Transcriptase Mutation Patterns..." (PMC2547475): positions
sufficient alone to influence at least one resistance-interpretation
algorithm. The "minor/accessory" list is deliberately NOT used as
primary ground truth here -- no single authoritative Stanford-sourced
list of exactly which positions count as minor could be confirmed from
available sources, and accessory mutations are frequently distal/
compensatory rather than direct pocket contacts, which this contact-
mining approach has no mechanism to detect regardless of accuracy.
"""

import os
import sys

import requests
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(__file__))
from mine_target import fetch_pdb, detect_ligand_resname
from prep_ensemble import prep_structure
from posegate.conserved_contacts import _structure_contact_residues, ensemble_reliability

UNIPROT_ACC = "P04585"
GAG_POL_OFFSET = 488  # verified constant: mature-protease author number + 488 = UniProt position
TOTAL_MONOMER_POSITIONS = 99

# RCSB-title-verified wild-type, ligand-bound HIV-1 protease structures.
# Excludes: explicit resistance mutants (1BV7/1BV9/1BWA/1BWB/1ODX/4Q1W-Y/
# 4Q5M), apo structures, and non-protease Gag-Pol domains (RNase H:
# 3QIN/3QIO; integrase: 5HRN/P/R/S, 5TC2) that share P04585's accession.
PDB_IDS = ["1DMP", "1HIV", "1HVH", "1HVR", "1HWR", "1ODY", "1QBR", "1QBS",
           "1QBT", "1QBU", "4U7Q", "5DGU", "5DGW", "5KAO"]

# Shafer, PMC2547475: positions sufficient alone to influence at least
# one PI resistance-interpretation algorithm rule.
MAJOR_RESISTANCE_POSITIONS = {30, 32, 46, 48, 50, 54, 82, 84, 88, 90}


def monomer_positions_contacted(structure) -> set:
    """This structure's specifically-contacted residues, collapsed from
    (chain, UniProt-remapped author number) onto a single monomer
    position 1-99, deduplicated so a residue contacted via BOTH chains
    in one structure counts once, not twice, toward that structure's
    contribution to the ensemble frequency."""
    residues = _structure_contact_residues(structure)
    if residues is None:
        return None
    positions = set()
    for r in residues:
        # r like 'GLU513.A' -> chain 'A', number 513
        chain = r[-1]
        num = int(''.join(c for c in r if c.isdigit()))
        positions.add(num - GAG_POL_OFFSET)
    return positions


def main():
    out_dir = "data/hiv_resistance_test"
    os.makedirs(out_dir, exist_ok=True)

    prepped = []
    for pdb_id in PDB_IDS:
        raw_path = fetch_pdb(pdb_id, out_dir)
        ligand_resname = detect_ligand_resname(raw_path)
        if ligand_resname is None:
            print(f"{pdb_id}: no ligand detected, skipping")
            continue
        try:
            s = prep_structure(pdb_id, raw_path, ligand_resname, out_dir, uniprot_acc=UNIPROT_ACC)
            prepped.append(s)
            print(f"{pdb_id}: prepped (ligand {ligand_resname})")
        except Exception as e:
            print(f"{pdb_id}: FAILED to prep ({e})")

    print(f"\nPrepped {len(prepped)}/{len(PDB_IDS)} structures.")

    # Per-structure monomer-collapsed contact sets, computed once and
    # reused for both the full-ensemble frequency and every LOO fold.
    per_structure_positions = {}
    for s in prepped:
        pos = monomer_positions_contacted(s)
        if pos is not None:
            per_structure_positions[s['pdb_id']] = pos
    usable = list(per_structure_positions.keys())
    n = len(usable)
    print(f"Usable (loaded + interaction detection succeeded): {n}/{len(prepped)}")

    # Consensus contact frequency per monomer position, across the full
    # (usable) ensemble.
    freq = {}
    for i in range(1, TOTAL_MONOMER_POSITIONS + 1):
        count = sum(1 for pid in usable if i in per_structure_positions[pid])
        if count > 0:
            freq[i] = count / n
    ranked = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))

    print(f"\nTop 20 monomer positions by consensus contact frequency (n_ensemble={n}):")
    print(f"{'pos':<6}{'freq':<8}{'n_structs':<10}resistance_position?")
    for i, f in ranked[:20]:
        flag = "YES" if i in MAJOR_RESISTANCE_POSITIONS else ""
        print(f"{i:<6}{f:<8.2f}{round(f*n):<10}{flag}")

    reliability = ensemble_reliability(n)
    print(f"\nEnsemble reliability: {reliability['tier'].upper()} -- {reliability['note']}")

    # LOO stability: for each held-out structure, does the SAME top-K set
    # (mined from the other n-1) actually appear as top-K again? Reported
    # as a stability fraction per K, not per residue -- a lighter-weight
    # signal than a full per-residue score, but grounded in the same
    # holdout logic as the rest of this project's LOO validation.
    for k in (10, 15):
        top_k_full = set(i for i, _ in ranked[:k])
        stability_hits = 0
        for held_pid in usable:
            train_pids = [p for p in usable if p != held_pid]
            fold_freq = {}
            for i in range(1, TOTAL_MONOMER_POSITIONS + 1):
                c = sum(1 for pid in train_pids if i in per_structure_positions[pid])
                if c > 0:
                    fold_freq[i] = c
            fold_top_k = set(i for i, _ in sorted(fold_freq.items(), key=lambda kv: (-kv[1], kv[0]))[:k])
            # A hit if the held-out structure's own contact set actually
            # touches at least one of the fold's top-K predicted positions.
            if fold_top_k & per_structure_positions[held_pid]:
                stability_hits += 1
        print(f"LOO top-{k} hit rate (held-out structure's own contacts include "
              f"a top-{k} prediction from the other {n-1}): {stability_hits}/{n} "
              f"= {stability_hits/n:.2f}")

    # Fisher's exact test, K=10 and K=15.
    print(f"\n{'='*70}\nFisher's exact test: Top-K mined positions vs "
          f"Stanford HIVdb major PI resistance positions\n{'='*70}")
    for k in (10, 15):
        top_k = set(i for i, _ in ranked[:k])
        a = len(top_k & MAJOR_RESISTANCE_POSITIONS)
        b = k - a
        c = len(MAJOR_RESISTANCE_POSITIONS) - a
        d = TOTAL_MONOMER_POSITIONS - a - b - c

        odds_ratio, p_two = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        _, p_greater = fisher_exact([[a, b], [c, d]], alternative='greater')

        recall = a / len(MAJOR_RESISTANCE_POSITIONS)
        expected_rate = len(MAJOR_RESISTANCE_POSITIONS) / TOTAL_MONOMER_POSITIONS
        observed_rate = a / k
        ef = observed_rate / expected_rate if expected_rate > 0 else float('nan')

        print(f"\nK={k}")
        print(f"  Top-{k}: {sorted(top_k)}")
        print(f"  Contingency table: [[a={a}, b={b}], [c={c}, d={d}]]")
        print(f"  Odds ratio:          {odds_ratio:.3f}")
        print(f"  p-value (two-sided): {p_two:.4f}")
        print(f"  p-value (greater):   {p_greater:.4f}")
        print(f"  Recall (a/(a+c)):    {recall:.3f}  ({a}/{len(MAJOR_RESISTANCE_POSITIONS)})")
        print(f"  Enrichment factor:   {ef:.2f}")


if __name__ == "__main__":
    main()
