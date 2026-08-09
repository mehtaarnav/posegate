# posegate/scripts/ca_family_matrix.py
"""Family-wide selectivity matrix across every human carbonic anhydrase
isoform with enough ligand-bound structures to mine.

The pairwise/three-way comparisons this project has run so far
(ca_family_selectivity.py, cdk_family_selectivity.py) answer "does this
residue differ between these two proteins". That is something a
medicinal chemist can already do by eye given the structures, which is
part of why a naive single-structure baseline recovers ~65% of the
literature selectivity residues (see BASELINE_COMPARISON_RESULT.md).

This asks the question that does NOT reduce to eyeballing a pair: across
an ENTIRE protein family at once, at every position that is a conserved
ligand contact in any member, what is the amino acid in every other
member? Positions where the identity is constant family-wide are the
shared catalytic scaffold -- useless for selectivity, and actively
dangerous to target. Positions where it varies are the selectivity
handles, and the matrix shows which specific isoforms they separate.

Output is one row per conserved-contact position (in reference-isoform
numbering) and one column per isoform, with the residue identity in
every cell and a marker for whether that isoform mined the position as
a top-N contact itself. Positions are sorted so the most variable
(most selectivity-relevant) appear first.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from posegate.selectivity import (
    fetch_uniprot_sequence, top_residue_numbers, build_alignment_map,
)

REFERENCE = 'CA2'
TOP_N = 10

# Only isoforms with enough ligand-bound structures to mine at all.
# Verified against UniProt (names and organism checked, not assumed) and
# RCSB structure counts: CA5A/CA5B have 0 ligand-bound entries <=2.5A,
# CA6 has 1, CA14 has 2, CA3 has 4 -- all excluded as unminable rather
# than mined at a size the reliability curve says is meaningless.
ISOFORMS = {
    'CA1':  {'acc': 'P00915', 'mined': 'data/ca1_verified/mined_result.json'},
    'CA2':  {'acc': 'P00918', 'mined': 'data/ca_verified/mined_result.json'},
    'CA4':  {'acc': 'P22748', 'mined': 'data/ca4_verified/mined_result.json'},
    'CA7':  {'acc': 'P43166', 'mined': 'data/ca7_verified/mined_result.json'},
    'CA9':  {'acc': 'Q16790', 'mined': 'data/ca9_verified/mined_result.json'},
    'CA12': {'acc': 'O43570', 'mined': 'data/ca12_verified/mined_result.json'},
    'CA13': {'acc': 'Q8N1Q1', 'mined': 'data/ca13_verified/mined_result.json'},
}


def main():
    available = {name: info for name, info in ISOFORMS.items()
                 if os.path.exists(info['mined'])}
    missing = set(ISOFORMS) - set(available)
    if missing:
        print(f"Skipping (not yet mined): {sorted(missing)}\n")
    if REFERENCE not in available:
        raise SystemExit(f"Reference isoform {REFERENCE} has no mined result; cannot build matrix.")

    sequences = {n: fetch_uniprot_sequence(i['acc']) for n, i in available.items()}
    top_native = {n: top_residue_numbers(i['mined'], TOP_N) for n, i in available.items()}
    ref_seq = sequences[REFERENCE]

    # Alignment maps in both directions against the reference isoform.
    to_ref, from_ref = {}, {}
    for name in available:
        if name == REFERENCE:
            continue
        m = build_alignment_map(ref_seq, sequences[name])
        from_ref[name] = m
        to_ref[name] = {v: k for k, v in m.items()}

    # Every reference-coordinate position that is a top-N conserved
    # contact in at least one isoform.
    #
    # A mined "residue" can be a retained metal ion rather than an amino
    # acid -- in this family that is the catalytic zinc (numbered past
    # the end of the protein sequence, e.g. ZN262 in a 260-residue CA2),
    # and it is a real, correctly-mined contact: every clinical CA
    # inhibitor is a zinc-binding sulfonamide. It has no sequence
    # position and no cross-isoform alignment, so it cannot appear in a
    # residue-identity matrix; it is collected and reported separately
    # rather than silently dropped, since dropping the single most
    # pharmacologically important contact in the family without comment
    # would be worse than not building the matrix at all.
    non_protein = {}
    positions = set()
    mined_by = {}

    def note(name, ref_pos):
        positions.add(ref_pos)
        mined_by.setdefault(ref_pos, set()).add(name)

    for native_pos in top_native[REFERENCE]:
        if native_pos > len(ref_seq):
            non_protein.setdefault(REFERENCE, []).append(native_pos)
        else:
            note(REFERENCE, native_pos)

    for name in available:
        if name == REFERENCE:
            continue
        for native_pos in top_native[name]:
            if native_pos > len(sequences[name]):
                non_protein.setdefault(name, []).append(native_pos)
                continue
            ref_pos = to_ref[name].get(native_pos)
            if ref_pos is None:
                continue
            note(name, ref_pos)

    names = sorted(available)
    rows = []
    for ref_pos in sorted(positions):
        cells = {}
        for name in names:
            if name == REFERENCE:
                aa = ref_seq[ref_pos - 1]
                native = ref_pos
            else:
                native = from_ref[name].get(ref_pos)
                aa = sequences[name][native - 1] if native else None
            cells[name] = (aa, native)
        identities = {aa for aa, _ in cells.values() if aa}
        rows.append({'ref_pos': ref_pos, 'cells': cells,
                      'n_identities': len(identities),
                      'mined_by': mined_by.get(ref_pos, set())})

    # Most variable positions first: those are the selectivity handles.
    rows.sort(key=lambda r: (-r['n_identities'], r['ref_pos']))

    print(f"Family-wide conserved-contact matrix, {len(names)} isoforms, "
          f"top-{TOP_N} each, positions in {REFERENCE} numbering")
    print(f"Cell = amino acid at the aligned position; * = that isoform mined "
          f"it as a top-{TOP_N} contact; '-' = no aligned residue\n")

    header = f"{'pos':>5}  {'var':>3}  " + "".join(f"{n:>7}" for n in names)
    print(header)
    print('-' * len(header))
    for r in rows:
        line = f"{r['ref_pos']:>5}  {r['n_identities']:>3}  "
        for name in names:
            aa, native = r['cells'][name]
            mark = '*' if name in r['mined_by'] else ' '
            line += f"{(aa or '-') + mark:>7}"
        print(line)

    print('-' * len(header))
    variable = [r for r in rows if r['n_identities'] > 1]
    invariant = [r for r in rows if r['n_identities'] == 1]
    print(f"\n{len(invariant)} invariant position(s) -- shared catalytic scaffold, "
          f"not selectivity-exploitable: {[r['ref_pos'] for r in invariant]}")
    print(f"{len(variable)} variable position(s) -- candidate selectivity handles, "
          f"most variable first: {[r['ref_pos'] for r in variable]}")

    if non_protein:
        print(f"\nNon-protein top-{TOP_N} contacts (retained metal ions -- real "
              f"contacts, but no sequence position, so not alignable into the "
              f"matrix above):")
        for name, nums in sorted(non_protein.items()):
            print(f"  {name}: {nums}")
        print("  For this family that is the catalytic zinc, which every clinical "
              "CA inhibitor coordinates -- pan-isoform by definition and therefore "
              "the opposite of a selectivity handle, but its absence from the "
              "matrix is a limitation of the representation, not evidence it "
              "does not matter.")


if __name__ == "__main__":
    main()
