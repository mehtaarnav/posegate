# posegate/scripts/cdk_family_selectivity.py
"""Cross-paralog selectivity comparison for CDK2 (P24941) vs CDK9
(P50750), the second family this method has been run on -- CA II/IX/XII
was the first (see ca_family_selectivity.py). The point of a second,
structurally unrelated family (kinase fold, not zinc-hydrolase) is to
check whether the CA result was a real capability or a fluke on one
enzyme class.

CDK2 vs CDK9 selectivity is a real, actively pursued oncology problem:
CDK9 (transcriptional, via P-TEFb) inhibitors are being developed for
cancer, and off-target CDK2 (cell-cycle) inhibition is a real selectivity
liability researchers explicitly try to avoid, and vice versa.

Residues are mapped between the two full-length UniProt sequences via
pairwise alignment (Bio.Align, BLOSUM62) rather than compared by raw
number, for the same reason as the CA family script: CDK2 and CDK9 are
different proteins with independent numbering, not two depositions of
the same protein.
"""

import json

import requests
from Bio import Align
from Bio.Align import substitution_matrices

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"

ISOFORMS = {
    'CDK2': {'acc': 'P24941', 'mined': 'data/cdk2_verified/mined_result.json'},
    'CDK9': {'acc': 'P50750', 'mined': 'data/cdk9_verified/mined_result.json'},
}
REFERENCE = 'CDK2'
TOP_N = 10


def fetch_sequence(acc: str) -> str:
    resp = requests.get(UNIPROT_FASTA_URL.format(acc=acc), timeout=30)
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    return ''.join(lines[1:])


def top_residue_numbers(mined_path: str, top_n: int) -> list:
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
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(ref_seq, other_seq)[0]

    mapping = {}
    for ref_block, other_block in zip(*alignment.aligned):
        for r, o in zip(range(*ref_block), range(*other_block)):
            mapping[r + 1] = o + 1
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
    other_name = [n for n in ISOFORMS if n != REFERENCE][0]
    other_seq = sequences[other_name]
    m_to_other = build_alignment_map(ref_seq, other_seq)
    m_to_ref = {v: k for k, v in m_to_other.items()}

    print(f"\nHinge-residue check (CDK2's top hit vs CDK9's top hit, mapped through the alignment):")
    cdk2_top1 = top['CDK2'][0]
    cdk9_top1 = top['CDK9'][0]
    mapped_from_cdk2 = m_to_other.get(cdk2_top1)
    mapped_from_cdk9 = m_to_ref.get(cdk9_top1)
    print(f"  CDK2 top-1 residue {cdk2_top1} ({ref_seq[cdk2_top1-1]}) aligns to "
          f"CDK9 position {mapped_from_cdk2} ({other_seq[mapped_from_cdk2-1] if mapped_from_cdk2 else '?'})")
    print(f"  CDK9 top-1 residue {cdk9_top1} ({other_seq[cdk9_top1-1]}) aligns to "
          f"CDK2 position {mapped_from_cdk9} ({ref_seq[mapped_from_cdk9-1] if mapped_from_cdk9 else '?'})")

    in_ref_coords = {REFERENCE: set(top[REFERENCE])}
    unmapped_report = {}
    mapped = set()
    unmapped = []
    for pos in top[other_name]:
        if pos in m_to_ref:
            mapped.add(m_to_ref[pos])
        else:
            unmapped.append(pos)
    in_ref_coords[other_name] = mapped
    if unmapped:
        print(f"\n{other_name}: {len(unmapped)}/{len(top[other_name])} top residue(s) unmapped "
              f"(gapped region): {unmapped}")

    print(f"\nTop-{TOP_N} residues, both expressed in {REFERENCE} reference numbering:")
    for name, positions in in_ref_coords.items():
        print(f"  {name}: {sorted(positions)}")

    shared = in_ref_coords[REFERENCE] & in_ref_coords[other_name]
    unique_ref = in_ref_coords[REFERENCE] - in_ref_coords[other_name]
    unique_other = in_ref_coords[other_name] - in_ref_coords[REFERENCE]
    print(f"\nShared (conserved kinase-domain scaffold candidates): {sorted(shared)}")
    print(f"Unique to {REFERENCE} (selectivity-relevant candidates): {sorted(unique_ref)}")
    print(f"Unique to {other_name} (selectivity-relevant candidates): {sorted(unique_other)}")

    print(f"\nNative residue identity at each shared/unique position (for literature cross-check):")
    for pos in sorted(shared | unique_ref | unique_other):
        ref_aa = ref_seq[pos-1]
        other_pos = m_to_other.get(pos)
        other_aa = other_seq[other_pos-1] if other_pos else '-'
        tag = ("shared" if pos in shared else
               f"{REFERENCE}-only" if pos in unique_ref else f"{other_name}-only")
        print(f"  CDK2-ref {pos}: CDK2={ref_aa}  CDK9={other_aa} ({'-' if other_pos is None else other_pos})  [{tag}]")


if __name__ == "__main__":
    main()
