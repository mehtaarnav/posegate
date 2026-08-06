# posegate/posegate/selectivity.py
"""Cross-paralog selectivity comparison: given several isoforms of a
protein family, each already mined independently by
scripts/mine_target.py (its own ligand-bound PDB ensemble, its own
UniProt accession, its own self-validated conserved-contact list), find
which top conserved-contact residues are shared across the family (the
scaffold) versus specific to one or a subset of isoforms (selectivity-
relevant candidates).

Extracted from two prior one-off scripts (CA II/IX/XII, then CDK2/CDK9 --
see conversation) into a single tested module once the same method had
independently reproduced a real, literature-confirmed selectivity
residue on two structurally unrelated families. Residue numbers are
compared through a pairwise sequence alignment against a chosen
reference isoform, never by raw number: different isoforms are different
proteins with independent numbering (e.g. CA9 carries an N-terminal PG
domain the others lack), and comparing raw numbers directly would repeat
the class of bug that caused the ERalpha 0%-accuracy failure (see
conversation), just between paralogs instead of between depositions of
one protein.
"""

from typing import Any, Dict, List, Optional

import requests
from Bio import Align
from Bio.Align import substitution_matrices

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def fetch_uniprot_sequence(acc: str) -> str:
    """The canonical sequence for a UniProt accession, as a plain string.
    Raises a clear error (not a bare requests exception) if the
    accession doesn't resolve or the response isn't FASTA -- an empty or
    malformed sequence here would otherwise fail confusingly deep inside
    the aligner instead of at the point of the actual problem."""
    try:
        resp = requests.get(UNIPROT_FASTA_URL.format(acc=acc), timeout=30)
        resp.raise_for_status()
    except requests.RequestException as e:
        raise ValueError(f"could not fetch UniProt sequence for accession {acc!r}: {e}") from e

    lines = resp.text.strip().splitlines()
    if not lines or not lines[0].startswith('>'):
        raise ValueError(f"UniProt response for {acc!r} was not FASTA "
                          f"(accession may not exist): {resp.text[:200]!r}")
    seq = ''.join(lines[1:])
    if not seq:
        raise ValueError(f"UniProt FASTA for {acc!r} had a header but no sequence")
    return seq


def top_residue_numbers(mined_path: str, top_n: int) -> List[int]:
    """Distinct residue numbers (int) from a mine_target.py output file,
    ranked by descending mined frequency, VdWContact excluded -- same
    convention as _top_k_predicted_residues in conserved_contacts.py.
    Returns [] (not an error) for a structure that mined nothing, since
    that's a legitimate, if uninformative, outcome for a small/failed
    ensemble -- callers that care should check len() themselves."""
    if top_n < 1:
        raise ValueError(f"top_n must be >= 1, got {top_n}")
    try:
        import json
        with open(mined_path) as f:
            data = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"{mined_path}: no mined_result.json here -- run scripts/mine_target.py "
            f"for this isoform first") from None
    except (json.JSONDecodeError, KeyError) as e:
        raise ValueError(f"{mined_path}: not a valid mine_target.py output ({e})") from e

    seen: List[int] = []
    for row in data.get('mined', []):
        if row['interaction'] == 'VdWContact':
            continue
        digits = ''.join(c for c in row['residue'] if c.isdigit())
        if not digits:
            continue  # malformed residue label; skip rather than crash on int()
        num = int(digits)
        if num not in seen:
            seen.append(num)
        if len(seen) >= top_n:
            break
    return seen


def build_alignment_map(ref_seq: str, other_seq: str) -> Dict[int, int]:
    """{ref_position (1-indexed): other_position (1-indexed)} for every
    alignment column with residues on both sides (match or mismatch, not
    a gap) -- a gapped position has no correspondence and is simply
    absent from the map, not guessed at."""
    if not ref_seq or not other_seq:
        raise ValueError("build_alignment_map requires two non-empty sequences")

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(ref_seq, other_seq)[0]

    mapping: Dict[int, int] = {}
    for ref_block, other_block in zip(*alignment.aligned):
        for r, o in zip(range(*ref_block), range(*other_block)):
            mapping[r + 1] = o + 1
    return mapping


def compare_isoforms(isoforms: Dict[str, Dict[str, str]], reference: str,
                      top_n: int = 10) -> Dict[str, Any]:
    """Compares top-N conserved-contact residues across N>=2 isoforms of
    a family, mapped onto `reference`'s numbering.

    Args:
        isoforms: {name: {'acc': UniProt accession, 'mined': path to
            that isoform's mine_target.py output}}, at least 2 entries.
        reference: which isoform's numbering to express every result in
            -- must be a key of `isoforms`.
        top_n: how many top mined residues per isoform to compare.

    Returns a dict with:
        'sequences': {name: sequence length}
        'top_native': {name: top-N residue numbers, that isoform's own numbering}
        'top_in_reference_coords': {name: set of positions mapped onto reference numbering}
        'unmapped': {name: [native positions with no reference-coordinate equivalent]}
        'shared_by_all': set of positions present in every isoform's top-N
        'pairwise_only': {(a, b): set of positions in exactly a's and b's top-N, no other}
        'unique_to': {name: set of positions only in that isoform's top-N}
    """
    if reference not in isoforms:
        raise ValueError(f"reference {reference!r} not in isoforms {list(isoforms)}")
    if len(isoforms) < 2:
        raise ValueError("compare_isoforms needs at least 2 isoforms to compare")

    sequences = {name: fetch_uniprot_sequence(info['acc']) for name, info in isoforms.items()}
    top_native = {name: top_residue_numbers(info['mined'], top_n) for name, info in isoforms.items()}

    ref_seq = sequences[reference]
    top_in_ref: Dict[str, set] = {reference: set(top_native[reference])}
    unmapped: Dict[str, list] = {}
    for name in isoforms:
        if name == reference:
            continue
        m_to_ref = {v: k for k, v in build_alignment_map(ref_seq, sequences[name]).items()}
        mapped, missed = set(), []
        for pos in top_native[name]:
            (mapped.add if pos in m_to_ref else missed.append)(m_to_ref[pos] if pos in m_to_ref else pos)
        top_in_ref[name] = mapped
        if missed:
            unmapped[name] = missed

    names = list(isoforms.keys())
    shared_by_all = set.intersection(*(top_in_ref[n] for n in names))

    pairwise_only = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            others = [n for n in names if n not in (a, b)]
            other_union = set().union(*(top_in_ref[n] for n in others)) if others else set()
            pairwise_only[(a, b)] = (top_in_ref[a] & top_in_ref[b]) - other_union

    unique_to = {}
    for name in names:
        other_union = set().union(*(top_in_ref[n] for n in names if n != name))
        unique_to[name] = top_in_ref[name] - other_union

    return {
        'sequences': {name: len(seq) for name, seq in sequences.items()},
        'top_native': top_native,
        'top_in_reference_coords': top_in_ref,
        'unmapped': unmapped,
        'shared_by_all': shared_by_all,
        'pairwise_only': pairwise_only,
        'unique_to': unique_to,
    }
