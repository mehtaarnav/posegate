# posegate/scripts/ca_family_selectivity.py
"""Cross-paralog selectivity comparison for the human carbonic anhydrase
family: CA II (P00918, the validated reference isoform), CA IX (Q16790,
tumor-associated), CA XII (O43570, tumor-associated).

CA II is the well-characterized, non-tumor "off-target" that most CA
inhibitors bind nonselectively; CA IX/XII selectivity over CA II is a
real, long-standing drug design problem. Each isoform is mined
independently (same pipeline, same self-validation, already proven on
5 unrelated targets -- see conversation) and the question here is purely
comparative: which conserved-contact residues are shared across all
three isoforms (the pan-CA catalytic scaffold), and which are specific
to one isoform (selectivity-relevant candidates)?

Residue numbers are NOT comparable across these three proteins directly
-- CA9 carries an N-terminal PG domain the others lack, so its catalytic
domain residues sit at a different absolute number even at the same
physical position. Comparing raw numbers here would repeat the same
mistake that caused the ERalpha 0%-accuracy bug (see conversation), just
between paralogs instead of between depositions of the same protein.
Residues are instead mapped through a pairwise sequence alignment of
each isoform's UniProt canonical sequence against CA II, using
Bio.Align.PairwiseAligner with BLOSUM62 -- standard for aligning
structurally near-identical, moderately-conserved catalytic domains.
"""

import json
import os

import requests
from Bio import Align
from Bio.Align import substitution_matrices

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

ISOFORMS = {
    'CA2':  {'acc': 'P00918', 'mined': 'data/ca_verified/mined_result.json'},
    'CA9':  {'acc': 'Q16790', 'mined': 'data/ca9_verified/mined_result.json'},
    'CA12': {'acc': 'O43570', 'mined': 'data/ca12_verified/mined_result.json'},
}
REFERENCE = 'CA2'
TOP_N = 10


def fetch_sequence(acc: str) -> str:
    resp = requests.get(UNIPROT_FASTA_URL.format(acc=acc), timeout=30)
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    return ''.join(lines[1:])  # drop the ">sp|..." header line


def top_residue_numbers(mined_path: str, top_n: int) -> list:
    """Distinct residue numbers (int) from a mine_target.py output,
    ranked by descending frequency, VdWContact excluded -- same
    convention as _top_k_predicted_residues in conserved_contacts.py."""
    with open(mined_path) as f:
        data = json.load(f)
    seen = []
    for row in data['mined']:
        if row['interaction'] == 'VdWContact':
            continue
        num = int(''.join(c for c in row['residue'] if c.isdigit()))
        if num not in seen:
            seen.append(num)
        if len(seen) >= top_n:
            break
    return seen


def build_alignment_map(ref_seq: str, other_seq: str) -> dict:
    """{ref_position (1-indexed) : other_position (1-indexed)} for every
    column the alignment calls a match/mismatch (not a gap on either
    side) -- a gapped position on either side has no correspondence and
    is simply absent from the map, not guessed."""
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(ref_seq, other_seq)[0]

    mapping = {}
    ref_pos, other_pos = 0, 0
    for ref_block, other_block in zip(*alignment.aligned):
        # aligned gives matched (start,end) index ranges (0-indexed) on
        # each sequence for each contiguous ungapped block.
        for r, o in zip(range(*ref_block), range(*other_block)):
            mapping[r + 1] = o + 1  # convert to 1-indexed UniProt numbering
    return mapping


def main():
    sequences = {name: fetch_sequence(info['acc']) for name, info in ISOFORMS.items()}
    for name, seq in sequences.items():
        print(f"{name} ({ISOFORMS[name]['acc']}): {len(seq)} residues")

    top = {name: top_residue_numbers(info['mined'], TOP_N) for name, info in ISOFORMS.items()}
    print(f"\nTop-{TOP_N} mined residues (own numbering, own protein):")
    for name, nums in top.items():
        print(f"  {name}: {nums}")

    ref_seq = sequences[REFERENCE]
    maps_from_ref = {}   # {other_name: {ref_pos: other_pos}}
    maps_to_ref = {}      # {other_name: {other_pos: ref_pos}}
    for name in ISOFORMS:
        if name == REFERENCE:
            continue
        m = build_alignment_map(ref_seq, sequences[name])
        maps_from_ref[name] = m
        maps_to_ref[name] = {v: k for k, v in m.items()}

    # Every isoform's top residues, expressed in CA2/reference numbering
    # where an alignment correspondence exists.
    in_ref_coords = {REFERENCE: set(top[REFERENCE])}
    for name in ISOFORMS:
        if name == REFERENCE:
            continue
        mapped = set()
        unmapped = []
        for pos in top[name]:
            if pos in maps_to_ref[name]:
                mapped.add(maps_to_ref[name][pos])
            else:
                unmapped.append(pos)
        in_ref_coords[name] = mapped
        if unmapped:
            print(f"\n{name}: {len(unmapped)}/{len(top[name])} top residue(s) fell in a "
                  f"gapped/unaligned region and have no CA2-coordinate equivalent: {unmapped}")

    print(f"\nTop-{TOP_N} residues, all expressed in {REFERENCE} reference numbering:")
    for name, positions in in_ref_coords.items():
        print(f"  {name}: {sorted(positions)}")

    shared_all_three = in_ref_coords['CA2'] & in_ref_coords['CA9'] & in_ref_coords['CA12']
    print(f"\nShared across all three isoforms (pan-CA catalytic scaffold candidates): "
          f"{sorted(shared_all_three)}")

    # Pairwise-only categories: a residue in exactly two isoforms' top-N
    # is neither "shared by all three" nor "unique to one" -- omitting
    # this category understated a real finding during development (see
    # conversation: CA2-ref position 130, the literature-documented CA
    # II/CA IX selectivity residue, sits here, shared by CA2+CA9 and
    # absent from CA12, which is exactly what the published mechanism
    # would predict -- it was invisible until this category existed).
    names = list(ISOFORMS.keys())
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            other = [n for n in names if n not in (a, b)][0]
            pair_only = (in_ref_coords[a] & in_ref_coords[b]) - in_ref_coords[other]
            print(f"Shared by {a}+{b} only (not {other}): {sorted(pair_only)}")

    for name in ISOFORMS:
        others = set().union(*(v for k, v in in_ref_coords.items() if k != name))
        unique = in_ref_coords[name] - others
        print(f"Unique to {name} (selectivity-relevant candidates): {sorted(unique)}")


if __name__ == "__main__":
    main()
