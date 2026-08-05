# posegate/scripts/resistance_position_test_pocket_baseline.py
"""Re-runs the resistance-position enrichment test (see
resistance_position_test.py) against a restricted denominator: the
active-site pocket (~20-30 residues) instead of the whole 99-residue
monomer.

Why: the original test's N=99 baseline can't distinguish "the miner
finds resistance-relevant residues specifically" from the much weaker
claim "the miner finds pocket residues, and pocket residues happen to be
where resistance mutations cluster" -- a real confound, since both the
mined contacts and the resistance positions are, by definition, near the
ligand. Restricting to a pocket-only denominator asks the harder
question: among residues already known to line the pocket, are the
miner's specific top picks still resistance-enriched, or is that where
the earlier signal actually lived?

The pocket must be defined independently of both things being compared,
or the test is circular:
  - NOT from the miner's own mined contacts (that's what's being tested)
  - NOT from the resistance position list (that's the ground truth)
  - NOT from literature recall (unverifiable, and this session already
    hit real errors trusting recalled PDB IDs and a recalled resistance
    list -- see conversation)

So the pocket is defined geometrically, by raw interatomic distance (not
ProLIF, to avoid sharing machinery with the method under test), from a
single reference structure -- 6DIF (RCSB-title-verified wild-type HIV-1
protease + tipranavir, same P04585-isolate accession as the mining
ensemble but not one of its 14 structures). Any monomer position with an
atom within 4.5A of any ligand atom, in either chain, counts as a
pocket residue.
"""

import os
import sys

import numpy as np
from rdkit import Chem
from scipy.stats import fisher_exact

sys.path.insert(0, os.path.dirname(__file__))
from mine_target import fetch_pdb, detect_ligand_resname
from prep_ensemble import prep_structure
from resistance_position_test import (
    PDB_IDS, UNIPROT_ACC, GAG_POL_OFFSET, MAJOR_RESISTANCE_POSITIONS,
    monomer_positions_contacted, TOTAL_MONOMER_POSITIONS,
)
from posegate.conserved_contacts import _structure_contact_residues, _load_structure

REFERENCE_PDB_ID = "6DIF"  # RCSB-title-verified wild-type protease + tipranavir,
                            # under the same P04585 accession as the mining
                            # ensemble, NOT one of its 14 structures
POCKET_CUTOFF_A = 4.5


def geometric_pocket_positions(structure) -> set:
    """Monomer positions (1-99, symmetry-collapsed) with any receptor
    atom within POCKET_CUTOFF_A of any ligand atom, computed by raw
    3D distance -- not ProLIF -- so this pocket definition shares no
    machinery with the interaction-detection method being tested."""
    lig, rec = _load_structure(structure)
    if lig is None or rec is None:
        raise ValueError(f"could not load reference structure {structure['pdb_id']}")

    lig_conf = lig.GetConformer()
    lig_pos = np.array([lig_conf.GetAtomPosition(i) for i in range(lig.GetNumAtoms())])

    rec_conf = rec.GetConformer()
    positions = set()
    for atom in rec.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is None:
            continue
        atom_pos = np.array(rec_conf.GetAtomPosition(atom.GetIdx()))
        if np.min(np.linalg.norm(lig_pos - atom_pos, axis=1)) <= POCKET_CUTOFF_A:
            author_num = pdb_info.GetResidueNumber()
            positions.add(author_num - GAG_POL_OFFSET)
    return positions


def main():
    out_dir = "data/hiv_resistance_test"
    os.makedirs(out_dir, exist_ok=True)

    print(f"Defining pocket from reference structure {REFERENCE_PDB_ID} "
          f"(not in the {len(PDB_IDS)}-structure mining ensemble)...")
    raw_path = fetch_pdb(REFERENCE_PDB_ID, out_dir)
    ligand_resname = detect_ligand_resname(raw_path)
    ref_structure = prep_structure(REFERENCE_PDB_ID, raw_path, ligand_resname, out_dir,
                                    uniprot_acc=UNIPROT_ACC)
    pocket = geometric_pocket_positions(ref_structure)
    print(f"Pocket (<= {POCKET_CUTOFF_A} A of ligand {ligand_resname} in "
          f"{REFERENCE_PDB_ID}): {len(pocket)} monomer positions -> {sorted(pocket)}")

    pocket_resistance_positions = pocket & MAJOR_RESISTANCE_POSITIONS
    print(f"Major resistance positions inside this pocket definition: "
          f"{sorted(pocket_resistance_positions)} ({len(pocket_resistance_positions)}/"
          f"{len(MAJOR_RESISTANCE_POSITIONS)} of the full resistance list)")

    # Re-mine the ensemble exactly as before (14 structures, all outside
    # the pocket-definition reference), to get the same ranked top-K.
    prepped = []
    for pdb_id in PDB_IDS:
        p_raw = fetch_pdb(pdb_id, out_dir)
        p_ligand = detect_ligand_resname(p_raw)
        if p_ligand is None:
            continue
        try:
            prepped.append(prep_structure(pdb_id, p_raw, p_ligand, out_dir, uniprot_acc=UNIPROT_ACC))
        except Exception as e:
            print(f"{pdb_id}: FAILED to prep ({e})")

    per_structure_positions = {}
    for s in prepped:
        pos = monomer_positions_contacted(s)
        if pos is not None:
            per_structure_positions[s['pdb_id']] = pos
    usable = list(per_structure_positions.keys())
    n = len(usable)

    freq = {}
    for i in range(1, TOTAL_MONOMER_POSITIONS + 1):
        count = sum(1 for pid in usable if i in per_structure_positions[pid])
        if count > 0:
            freq[i] = count / n
    ranked = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0]))

    print(f"\n{'='*70}\nFisher's exact test, restricted to the {len(pocket)}-residue "
          f"geometric pocket (not N={TOTAL_MONOMER_POSITIONS})\n{'='*70}")
    for k in (10, 15):
        top_k_all = set(i for i, _ in ranked[:k])
        # Restricted to the pocket universe: positions outside the
        # pocket definition are excluded entirely, not counted as "not
        # in top-K" against the full monomer -- that's the actual fix.
        top_k_in_pocket = top_k_all & pocket

        a = len(top_k_in_pocket & MAJOR_RESISTANCE_POSITIONS)
        b = len(top_k_in_pocket) - a
        c = len(pocket_resistance_positions) - a
        d = len(pocket) - a - b - c

        if a + b == 0 or c + d == 0:
            print(f"\nK={k}: degenerate table (top-K entirely outside pocket, or "
                  f"pocket has no non-resistance residues) -- skipping")
            continue

        odds_ratio, p_two = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        _, p_greater = fisher_exact([[a, b], [c, d]], alternative='greater')
        recall = a / len(pocket_resistance_positions) if pocket_resistance_positions else float('nan')
        expected_rate = len(pocket_resistance_positions) / len(pocket)
        observed_rate = a / len(top_k_in_pocket) if top_k_in_pocket else 0.0
        ef = observed_rate / expected_rate if expected_rate > 0 else float('nan')

        print(f"\nK={k} (of which {len(top_k_in_pocket)}/{k} fall inside the pocket "
              f"definition; {k - len(top_k_in_pocket)} mined position(s) outside it, "
              f"excluded from this test)")
        print(f"  Top-K in pocket:     {sorted(top_k_in_pocket)}")
        print(f"  Pocket universe:     {len(pocket)} residues, of which "
              f"{len(pocket_resistance_positions)} are major resistance positions")
        print(f"  Contingency table:   [[a={a}, b={b}], [c={c}, d={d}]]")
        print(f"  Odds ratio:          {odds_ratio:.3f}")
        print(f"  p-value (two-sided): {p_two:.4f}")
        print(f"  p-value (greater):   {p_greater:.4f}")
        print(f"  Recall (a/(a+c)):    {recall:.3f}")
        print(f"  Enrichment factor:   {ef:.2f}  (vs pocket-internal baseline rate "
              f"{expected_rate:.2f}, not the whole-monomer rate)")


if __name__ == "__main__":
    main()
