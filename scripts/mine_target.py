# posegate/scripts/mine_target.py
"""End-to-end entry point: give it PDB IDs, get a mined, self-validated
conserved contact back. No manifest to hand-build, no ligand resname to
look up, no WSL-specific orchestration script.

    python scripts/mine_target.py --pdb_ids 4FKL 1GZ8 5JQ5 3EZR 3WBL \
        --out_dir data/my_target --uniprot_acc P00533

--uniprot_acc is required, not optional: found directly on ERalpha that a
keyword-based candidate ID list can silently pull in structures of a
different protein (ERbeta, ERR-gamma both slipped into an "ERalpha" list
this way), and that even genuine same-protein depositions can use
inconsistent author residue numbering across entries -- both of which
corrupt mining without ever raising an error. Every fetched structure is
checked against this accession via SIFTS and remapped onto consistent
UniProt numbering before mining; a structure that maps to a different
accession, or that SIFTS has no mapping for, is rejected rather than
silently pooled in wrong. See posegate.residue_mapping.

Everything else this project has done so far -- the five-family
validation, the visGReMLIN comparison, the LOO self-validation -- was
run through prep_ensemble.py against a hand-built {pdb_id, pdb_path,
ligand_resname} manifest, which requires knowing each structure's bound
ligand's three-letter PDB code in advance. That is exactly the kind of
domain knowledge a researcher new to a target does not yet have, so this
script infers it: for each fetched structure, it picks the largest
non-solvent, non-buffer HETATM residue as the ligand. This is a
heuristic, not a guarantee -- see LIKELY_NON_LIGAND below -- and is
reported per structure so a wrong guess is visible, not silent.
"""

import argparse
import json
import os
from datetime import datetime, timezone

import requests

from posegate.conserved_contacts import mine_conserved_contacts, leave_one_out_validate
from posegate.numbering_display import residue_numbering_table
# Sibling script, not a package module (scripts/ has no __init__.py):
# importable directly because running `python scripts/mine_target.py`
# puts this file's own directory on sys.path.
from prep_ensemble import prep_structure

# Common crystallization buffer components, cryoprotectants and metal
# ions that appear as HETATM records but are not the bound ligand a
# conserved-contact analysis cares about. Not exhaustive -- a structure
# whose real ligand is smaller than one of these, or an unlisted buffer
# component, will be misdetected, which is why the chosen ligand is
# printed for every structure rather than assumed correct silently.
LIKELY_NON_LIGAND = {
    'HOH', 'SO4', 'GOL', 'EDO', 'DMS', 'PEG', 'PO4', 'CL', 'NA', 'K', 'MG',
    'CA', 'ZN', 'ACT', 'TRS', 'BME', 'IMD', 'MPD', '1PE', 'PGE', 'MRD',
    'FMT', 'EPE', 'NO3', 'IOD', 'BR', 'FLC', 'CIT', 'TLA', 'MLI', 'ACY',
    'NH4', 'SCN', 'AZI', 'P6G', '2PE', 'PG4', '12P', '15P', 'XPE', 'DTT',
    'TCE', 'CO', 'MN', 'NI', 'CD', 'CS', 'RB', 'SR', 'BA', 'LI', 'F',
    'OXY', 'UNX', 'UNL', 'MN3', 'CU', 'CU1', 'FE', 'FE2',
    # Heavy-atom isomorphous-replacement/anomalous-phasing derivatives,
    # common in older crystal structures, retained in the deposited
    # coordinates but playing no role in binding. Found the hard way: a
    # bound mercury ion got picked as "the ligand" for two carbonic
    # anhydrase structures, and ProLIF's VdWContact then crashed outright
    # (no van der Waals radius for Hg in its chosen radii table) rather
    # than just producing a nonsensical result.
    'HG', 'PT', 'AU', 'AG', 'OS', 'IR', 'PD', 'W', 'RE', 'SM', 'GD', 'YB',
    'TB', 'EU', 'LU', 'PB', 'U', 'TH', 'HO',
    # Prosthetic-group cofactors: permanently bound, structural/catalytic
    # (e.g. heme's iron enables COX's peroxidase activity), not the
    # orthosteric small-molecule ligand a conserved-contact analysis
    # cares about -- and, being large, heavy-atom-count residues, they
    # beat real drug-like ligands under the "largest HETATM residue"
    # heuristic. Found the hard way: HEM (heme B) was picked as "the
    # ligand" for every COX-1 structure in this ensemble (COX is a
    # heme-dependent bifunctional peroxidase/cyclooxygenase), crashing
    # the whole run outright when RDKit/obabel couldn't cleanly convert
    # heme's porphyrin ring system through the ligand-prep pipeline.
    'HEM', 'HEC', 'HEA', 'HEB', 'HEG', 'DHE', 'HDD', 'FAD', 'FMN', 'NAD',
    'NAP', 'NDP', 'FES', 'SF4', 'F3S', 'BCL', 'CLA', 'PLP',
    # Heme-derivative variants (metal-substituted or metal-free
    # porphyrin) found the same way as HEM itself on COX-1, a heme-
    # dependent peroxidase: COH (metal-free protoporphyrin IX), MNH
    # (manganese-reconstituted protoporphyrin IX, used in some COX
    # structures to trap a distinct catalytic state).
    'COH', 'MNH', 'PP9', 'DDH',
    # Detergents used to solubilize membrane-associated proteins for
    # crystallization (COX is membrane-associated) -- structurally
    # present in the crystal but not a drug-like orthosteric ligand.
    'BOG', 'LDA', 'OGA', 'HTG', 'C8E', 'DDQ', 'LMT', 'NG6',
    # N/O-linked glycosylation sugars: present on essentially every
    # secreted or membrane glycoprotein (e.g. COX-2), not a drug ligand.
    # Found the same way on ERalpha earlier this session (see
    # conversation) and again here on COX-2, where NAG outcompeted the
    # real bound inhibitor/substrate by atom count in 3 of 7 structures.
    'NAG', 'MAN', 'BMA', 'FUC', 'GAL', 'NDG', 'NGA', 'SIA', 'XYP',
}


def fetch_pdb(pdb_id: str, out_dir: str) -> str:
    path = os.path.join(out_dir, f"{pdb_id}.pdb")
    if os.path.exists(path):
        return path
    resp = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=30)
    resp.raise_for_status()
    with open(path, 'w') as f:
        f.write(resp.text)
    return path


# Above this many PDB entries containing a given chemical component
# anywhere, it's treated as a crystallization additive/cofactor/glycan
# rather than a specific bound ligand, regardless of atom count. Set
# from a real observed sample (see global_pdb_prevalence's docstring):
# known additives ranged 515-28155 entries, known real ligands 6-94 --
# 200 sits cleanly between them with margin on both sides.
PREVALENCE_THRESHOLD = 200
_prevalence_cache: dict = {}


def global_pdb_prevalence(comp_id: str, timeout: int = 15) -> int:
    """How many PDB entries contain this chemical component anywhere,
    queried live from RCSB. A drug-like ligand specific to one binding
    site appears in a handful of entries; a crystallization additive,
    cofactor, or glycosylation sugar appears in hundreds to tens of
    thousands, essentially regardless of the protein. This is the
    general, self-updating replacement for LIKELY_NON_LIGAND's manually-
    curated list, which cannot cover every contaminant class in advance
    -- three different classes (heme cofactors, detergents, glycans)
    each needed their own hand-added entries this session before this
    function existed, and there is no reason a fourth, unseen class
    wouldn't need one too, indefinitely, without this fix.

    Verified on a real sample: HEM 6493, NAG 12288, SO4 28155, GOL 26103,
    BOG 515 entries, versus real bound ligands CEL 6, STU 94, IBP 14 --
    a wide, clean separation (see conversation).

    Cached per process (an ensemble mines the same handful of comp_ids
    repeatedly). Returns -1, not an exception, if the query fails --
    callers should treat that as 'unknown' and fail open (allow the
    candidate) rather than block a whole run over a network hiccup.
    """
    if comp_id in _prevalence_cache:
        return _prevalence_cache[comp_id]
    try:
        query = {
            'query': {'type': 'terminal', 'service': 'text', 'parameters': {
                'attribute': 'rcsb_nonpolymer_instance_annotation.comp_id',
                'operator': 'exact_match', 'value': comp_id}},
            'return_type': 'entry',
            'request_options': {'return_counts': True},
        }
        resp = requests.post('https://search.rcsb.org/rcsbsearch/v2/query',
                              json=query, timeout=timeout)
        if resp.status_code == 204:
            count = 0
        else:
            resp.raise_for_status()
            count = resp.json().get('total_count', -1)
    except requests.RequestException:
        count = -1
    _prevalence_cache[comp_id] = count
    return count


def detect_ligand_resname(pdb_path: str):
    """The largest (by atom count) HETATM residue instance in the file
    that is neither a known non-ligand (LIKELY_NON_LIGAND, a fast
    offline first pass for the most common, already-catalogued cases)
    nor globally common across the whole PDB (global_pdb_prevalence,
    checked live -- see its docstring for why this generalizes beyond
    what LIKELY_NON_LIGAND can cover). Returns None if nothing qualifies.
    """
    counts = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            resname = line[17:20].strip()
            if resname in LIKELY_NON_LIGAND:
                continue
            key = (resname, line[21], line[22:26])
            counts[key] = counts.get(key, 0) + 1
    if not counts:
        return None

    # Largest single residue instance per resname, not summed across
    # chains/copies -- a resname repeated across many small
    # crystallographic copies should not outrank one real, larger,
    # single bound ligand.
    best_per_resname: dict = {}
    for (resname, _, _), n in counts.items():
        best_per_resname[resname] = max(best_per_resname.get(resname, 0), n)

    for resname, _ in sorted(best_per_resname.items(), key=lambda kv: -kv[1]):
        prevalence = global_pdb_prevalence(resname)
        # prevalence == -1 means the query failed (network issue, not a
        # verdict) -- fail open rather than block the whole run.
        if prevalence == -1 or prevalence <= PREVALENCE_THRESHOLD:
            return resname
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_ids", nargs='+', required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--top_n", type=int, default=10)
    parser.add_argument("--skip_self_validation", action="store_true")
    parser.add_argument("--uniprot_acc", required=True,
                         help="UniProt accession (e.g. P03372) every --pdb_ids entry is expected "
                              "to be. Required, not optional: found directly on ERalpha that a "
                              "keyword-based candidate ID list can silently include structures of "
                              "a different protein entirely (ERbeta, ERR-gamma both slipped into "
                              "an 'ERalpha' list), and that even genuine ERalpha depositions use "
                              "inconsistent author residue numbering across entries -- both of "
                              "which silently corrupt mining without ever raising an error. Every "
                              "fetched structure is checked against this accession via SIFTS and "
                              "remapped onto its numbering before mining; a structure that maps to "
                              "a different accession is rejected, not silently pooled in wrong.")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    prepped = []
    for pdb_id in args.pdb_ids:
        print(f"{pdb_id}: fetching...")
        try:
            raw_path = fetch_pdb(pdb_id, args.out_dir)
        except requests.RequestException as e:
            print(f"  FAILED to fetch: {e}")
            continue

        ligand_resname = detect_ligand_resname(raw_path)
        if ligand_resname is None:
            print(f"  FAILED: no non-solvent/buffer HETATM residue found "
                  f"(apo structure, or its ligand is in LIKELY_NON_LIGAND -- check manually)")
            continue
        print(f"  detected ligand: {ligand_resname} (verify this is right -- see the module "
              f"docstring on how this is guessed)")

        try:
            prepped.append(prep_structure(pdb_id, raw_path, ligand_resname, args.out_dir,
                                           uniprot_acc=args.uniprot_acc))
        except Exception as e:
            print(f"  FAILED to prep: {e}")

    print(f"\nPrepped {len(prepped)}/{len(args.pdb_ids)} structures.")
    if len(prepped) < 2:
        raise SystemExit("Need at least 2 successfully prepped structures to mine anything.")

    asymmetric_ids = [s['pdb_id'] for s in prepped if s.get('asymmetric_multichain')]
    if asymmetric_ids:
        print(f"\n{'!'*70}")
        print(f"UNSUPPORTED TARGET CLASS WARNING: {len(asymmetric_ids)}/{len(prepped)} structure(s) "
              f"have a receptor with multiple, non-identical chains near the ligand: {asymmetric_ids}")
        print("This pipeline has no mechanism to verify chain-letter assignment is consistent "
              "for a non-symmetric multi-chain protein across different PDB depositions -- found "
              "on chymotrypsin (three fragments of one cleaved polypeptide), where inconsistent "
              "chain labeling fragmented the same physical residue's identity across the ensemble "
              "and produced an unusable signal (0% self-validation accuracy) despite looking like "
              "a normal run. A true symmetric multimer (e.g. HIV protease's homodimer, both chains "
              "identical) does not have this problem and is not flagged here.")
        print("Results below should be treated as UNVERIFIED for this reason, not just low-confidence.")
        print(f"{'!'*70}")

    print(f"\nMining conserved contacts across {len(prepped)} structures...")
    results = mine_conserved_contacts(prepped)
    specific = [r for r in results if r['interaction'] != 'VdWContact']
    top_specific = residue_numbering_table(specific, prepped, top_n=args.top_n)

    # 'Native #' is the author residue number(s) this row's UniProt
    # position actually corresponds to in the source structures -- what
    # a user needs to type into PyMOL/ChimeraX, not the UniProt number
    # in the 'Residue' column. Shown side by side so a raw-number
    # assumption never has to be checked with a one-off script again
    # (see posegate.numbering_display's docstring for why this exists).
    print(f"\n{'Residue':<14}{'Interaction':<14}{'N':>4}{'Frequency':>12}  95% CI     {'Native #'}")
    print('-' * 80)
    for r in top_specific:
        lo, hi = r['ci95']
        print(f"{r['residue']:<14}{r['interaction']:<14}{r['n_structures']:>4}{r['frequency']:>12.2f}"
              f"  [{lo:.2f}, {hi:.2f}]  {r['native_numbering']}")

    self_validation = None
    if not args.skip_self_validation:
        print(f"\nSelf-validating: leave-one-out over these same {len(prepped)} structures...")
        self_validation = leave_one_out_validate(prepped)
        print(f"\n{'top-k':<8}{'hits/folds':<14}accuracy")
        print('-' * 32)
        for k, acc in self_validation['accuracy'].items():
            val = f"{acc['accuracy']:.2f}" if acc['accuracy'] is not None else "n/a"
            hits_str = f"{acc['hits']}/{acc['n_folds']}"
            print(f"top-{k:<4}{hits_str:<14}{val}")
        reliability = self_validation['reliability']
        print(f"\nEnsemble size reliability: {reliability['tier'].upper()}")
        print(f"  {reliability['note']}")

    provenance = {
        'requested_pdb_ids': args.pdb_ids,
        'prepped_pdb_ids': [s['pdb_id'] for s in prepped],
        'uniprot_acc': args.uniprot_acc,
        'top_n': args.top_n,
        'generated_at_utc': datetime.now(timezone.utc).isoformat(),
        'asymmetric_multichain_structures': asymmetric_ids,
        'unsupported_target_class': bool(asymmetric_ids),
    }

    # Every mined row, not just the printed top-N, carries its native
    # (author) numbering in the saved JSON -- a downstream reader
    # shouldn't need to recompute this from the prepped structures to
    # answer "what do I actually select in PyMOL for this residue."
    mined_with_native = residue_numbering_table(results, prepped, top_n=None)

    out_json = os.path.join(args.out_dir, 'mined_result.json')
    with open(out_json, 'w') as f:
        json.dump({'provenance': provenance, 'mined': mined_with_native,
                    'self_validation': self_validation}, f, indent=2)
    print(f"\nWrote {out_json}")


if __name__ == "__main__":
    main()
