# posegate/posegate/numbering_display.py
"""Shows a mined residue's UniProt numbering (the internal key
everything in this project uses) alongside the actual author numbering
it corresponds to in each source structure -- the number a user would
need to type into PyMOL/ChimeraX to select that residue, which is
usually NOT the same as the UniProt number displayed everywhere else.

Exists because this project has repeatedly needed a one-off diagnostic
script to answer "what does this UniProt-numbered residue actually
correspond to in the original depositions" -- most consequentially, a
wrong assumption that raw SIFTS-remapped numbers tracked classical
literature numbering caused a real misdiagnosis on the trypsin hold-out
test (reported as "falsified" and later corrected -- see
HOLDOUT_RESULT_trypsin_chymotrypsin.md). If every mined result had shown
its author numbering alongside the UniProt number from the start, that
mistake would have been visibly avoidable rather than requiring a
separate motif-anchoring investigation after the fact.

A second, independent value: if a single UniProt position maps to
DIFFERENT author numbers across structures that are supposed to be
consistent, that inconsistency is directly visible here rather than
silently fragmenting the mined signal the way the ERalpha numbering bug
did before SIFTS remapping existed at all (see posegate.residue_mapping).
Consistency across the ensemble is not assumed; it's shown.
"""

from typing import Any, Dict, List


def build_multi_coordinate_index(structures: List[Dict[str, Any]]) -> Dict[int, Dict[str, str]]:
    """{uniprot_position: {pdb_id: 'chain+author_number'}} built from
    every structure's own residue_map (see prep_ensemble.prep_structure).
    Structures with no residue_map (uniprot_acc wasn't used) contribute
    nothing -- this display has nothing to show without it."""
    index: Dict[int, Dict[str, str]] = {}
    for s in structures:
        residue_map = s.get('residue_map')
        if not residue_map:
            continue
        for (chain, author_num), unp_pos in residue_map.items():
            index.setdefault(unp_pos, {})[s['pdb_id']] = f"{chain}{author_num}"
    return index


def format_native_numbering(index: Dict[int, Dict[str, str]], uniprot_pos: int) -> str:
    """Human-readable native/author numbering for one UniProt position
    across the ensemble. If every structure that covers this position
    agrees on the author number, shown once (e.g. '189'). If they
    disagree, every distinct variant is shown with a flag -- the
    disagreement itself is the finding, not something to average away.
    '?' if no structure in the index covers this position at all."""
    entries = index.get(uniprot_pos, {})
    if not entries:
        return "?"
    distinct = sorted(set(entries.values()))
    if len(distinct) == 1:
        return distinct[0]
    return "/".join(distinct) + " (INCONSISTENT across ensemble)"


def residue_numbering_table(mined_rows: List[Dict[str, Any]], structures: List[Dict[str, Any]],
                             top_n: int = None) -> List[Dict[str, Any]]:
    """Augments mined_rows (mine_conserved_contacts' output) with a
    'native_numbering' field per row, without mutating the input.
    top_n limits to the first N rows (mined_rows is expected pre-sorted
    by frequency, matching mine_conserved_contacts' contract) -- pass
    None for every row."""
    index = build_multi_coordinate_index(structures)
    rows = mined_rows[:top_n] if top_n is not None else mined_rows
    out = []
    for row in rows:
        digits = ''.join(c for c in row['residue'] if c.isdigit())
        unp_pos = int(digits) if digits else None
        native = format_native_numbering(index, unp_pos) if unp_pos is not None else "?"
        out.append({**row, 'native_numbering': native})
    return out
