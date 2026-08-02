# posegate/scripts/compare_visgremlin.py
"""Head-to-head of posegate's conserved-contact miner against visGReMLIN
[Ribeiro et al. 2020, BMC Bioinformatics 21:80] on visGReMLIN's own
published CDK case study.

visGReMLIN was released as a web server only. No source code or package
was published, and as of this writing both advertised URLs are
unreachable (vagner.dti.ufv.br/visgremlin4 refuses connections;
homepages.dcc.ufmg.br/~alexandrefassio/gremlin/ returns PHP errors), so
it cannot be run on new data. This script compares against its published
results on its own CDK case study, which scores motif recovery against
the experimentally determined CDK binding-site atoms of Schonbrunn et
al. Those atoms, and visGReMLIN's reported recovery of 18 of 26 (69%),
are transcribed from Table 2 of the paper below.

Two asymmetries follow from this setup and are reported alongside the
scores:

  1. Ensemble size. visGReMLIN mined 73 CDK-inhibitor complexes, while
     the ensembles here are 6 to 22 structures. A smaller ensemble is a
     harder setting for a frequency-based method, since each structure
     carries proportionally more weight in the frequency, but it is not
     an equivalent one.
  2. Granularity. visGReMLIN's motifs are sets of atoms, so its 18/26 is
     an atom-level score. posegate aggregates per residue and produces no
     atom-level output, so it cannot be scored on that denominator. The
     comparable score is residue-level: of the 9 reference residues, how
     many does each method place at the interface? visGReMLIN's
     residue-level score is derived from Table 2 by counting a residue as
     recovered if any of its atoms is marked, which is the most favorable
     reading available for it.
"""

import argparse
import json
import re
from collections import defaultdict

# visGReMLIN 2020, Table 2: the Schonbrunn et al. experimentally
# determined CDK binding-site atoms, and whether a visGReMLIN motif
# contained each. The paper marks cells with a check (found), a dot, or a
# cross; only checks are counted in its own stated 18/26, so both other
# markers are recorded as not-found here.
VISGREMLIN_TABLE2 = {
    'ASP145': {'CB': True, 'CG': False, 'OD1': True},
    'LYS33': {'CB': True, 'CD': True, 'CE': False, 'CG': True, 'NZ': False},
    'ASP86': {'N': True, 'CB': True, 'OD1': False, 'OD2': True},
    'LYS89': {'CB': True, 'CE': False, 'NZ': False},
    'HIS84': {'O': False},
    'LEU83': {'N': True, 'O': True},
    'PHE82': {'CE2': True, 'CZ': True},
    'GLU81': {'O': True},
    'PHE80': {'CB': True, 'CG': True, 'CD2': True, 'CE2': True, 'CZ': False},
}

# ProLIF's VdWContact fires on any pair of atoms inside their van der
# Waals radii, and is not an interaction type visGReMLIN models: its
# vocabulary is aromatic stacking, hydrogen bond, hydrophobic, repulsive
# and salt bridge. Including VdWContact therefore credits posegate for a
# criterion the comparator never reported, so both variants are printed
# and the specific-only variant is the conservative one.
POSEGATE_SPECIFIC = {'HBDonor', 'HBAcceptor', 'Hydrophobic', 'FaceToFace', 'EdgeToFace'}


def residue_key(residue: str) -> str:
    """ProLIF ResidueId repr ('LEU83.A') -> reference-set key ('LEU83')."""
    return re.match(r'^([A-Za-z]+\d+)', residue).group(1).upper()


def best_frequencies(freq_table: list, specific_only: bool) -> dict:
    """residue -> highest frequency at which the miner reports it, and the
    interaction type carrying that frequency."""
    out = defaultdict(lambda: (0.0, None))
    for row in freq_table:
        if specific_only and row['interaction'] not in POSEGATE_SPECIFIC:
            continue
        key = residue_key(row['residue'])
        if row['frequency'] > out[key][0]:
            out[key] = (row['frequency'], row['interaction'])
    return dict(out)


def score(freqs: dict, threshold: float) -> set:
    return {r for r in VISGREMLIN_TABLE2 if freqs.get(r, (0.0, None))[0] >= threshold}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq_json", required=True,
                        help="Miner output for a CDK2 ensemble (run_conserved_contact_miner.py --out_json)")
    parser.add_argument("--thresholds", nargs='+', type=float, default=[0.0001, 0.33, 0.5, 0.67],
                        help="Frequency thresholds at which to score posegate")
    parser.add_argument("--out_json", default=None)
    args = parser.parse_args()

    with open(args.freq_json) as f:
        table = json.load(f)

    vg_residues = {r for r, atoms in VISGREMLIN_TABLE2.items() if any(atoms.values())}
    vg_atoms_found = sum(sum(a.values()) for a in VISGREMLIN_TABLE2.values())
    vg_atoms_total = sum(len(a) for a in VISGREMLIN_TABLE2.values())
    n_ref = len(VISGREMLIN_TABLE2)

    print("Reference set: Schonbrunn et al. CDK binding site, as tabulated by "
          "visGReMLIN (Ribeiro et al. 2020, Table 2)")
    print(f"  {vg_atoms_total} atoms across {n_ref} residues\n")
    print(f"visGReMLIN, atom level:    {vg_atoms_found}/{vg_atoms_total} "
          f"({vg_atoms_found / vg_atoms_total:.0%})  [as published; 73-complex CDK ensemble]")
    print(f"visGReMLIN, residue level: {len(vg_residues)}/{n_ref} "
          f"({len(vg_residues) / n_ref:.0%})  [residue counted found if any atom checked]")
    print(f"  missed entirely: {sorted(set(VISGREMLIN_TABLE2) - vg_residues)}\n")

    results = {'visgremlin': {'atoms_found': vg_atoms_found, 'atoms_total': vg_atoms_total,
                              'residues_found': sorted(vg_residues), 'n_ref_residues': n_ref},
               'posegate': {}}

    for specific_only in (True, False):
        label = "specific interactions only" if specific_only else "including VdWContact"
        freqs = best_frequencies(table, specific_only)
        print(f"posegate, residue level ({label}):")
        for t in args.thresholds:
            found = score(freqs, t)
            shown = 'any' if t <= 0.0001 else f'>={t:.2f}'
            print(f"  freq {shown:<8} {len(found)}/{n_ref} ({len(found) / n_ref:.0%})"
                  f"   missed: {sorted(set(VISGREMLIN_TABLE2) - found)}")
            results['posegate'].setdefault(label, {})[shown] = sorted(found)
        print()

    print(f"Per-residue detail (posegate, best frequency over the ensemble):")
    print(f"  {'residue':<10}{'visGReMLIN':<12}{'posegate freq':>14}  interaction")
    all_freqs = best_frequencies(table, specific_only=False)
    spec_freqs = best_frequencies(table, specific_only=True)
    for r in sorted(VISGREMLIN_TABLE2, key=lambda x: -all_freqs.get(x, (0.0, None))[0]):
        f_any, i_any = all_freqs.get(r, (0.0, None))
        f_spec, i_spec = spec_freqs.get(r, (0.0, None))
        vg = 'found' if r in vg_residues else 'missed'
        detail = f"{i_any}" if i_any else '-'
        if i_spec and i_spec != i_any:
            detail += f" (specific: {i_spec} @ {f_spec:.2f})"
        elif not i_spec:
            detail += " (no specific contact)"
        print(f"  {r:<10}{vg:<12}{f_any:>14.2f}  {detail}")

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
