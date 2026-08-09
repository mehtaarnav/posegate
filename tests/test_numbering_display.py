# posegate/tests/test_numbering_display.py
from posegate.numbering_display import (
    build_multi_coordinate_index, format_native_numbering, residue_numbering_table,
)


def _structure(pdb_id, residue_map):
    return {'pdb_id': pdb_id, 'residue_map': residue_map}


def test_build_multi_coordinate_index_maps_unp_position_to_author_per_structure():
    structures = [
        _structure('1AAA', {('A', 300): 500, ('A', 301): 501}),
        _structure('1BBB', {('A', 1): 500, ('A', 2): 501}),
    ]
    index = build_multi_coordinate_index(structures)
    assert index[500] == {'1AAA': 'A300', '1BBB': 'A1'}
    assert index[501] == {'1AAA': 'A301', '1BBB': 'A2'}


def test_build_multi_coordinate_index_skips_structures_without_a_residue_map():
    structures = [
        _structure('1AAA', {('A', 300): 500}),
        {'pdb_id': '1CCC', 'residue_map': None},
        {'pdb_id': '1DDD'},  # key entirely absent
    ]
    index = build_multi_coordinate_index(structures)
    assert index == {500: {'1AAA': 'A300'}}


def test_format_native_numbering_consistent_case():
    index = {500: {'1AAA': 'A194', '1BBB': 'A194', '1CCC': 'A194'}}
    assert format_native_numbering(index, 500) == 'A194'


def test_format_native_numbering_flags_inconsistency_across_ensemble():
    """The real motivating case: the same UniProt position corresponds
    to different author numbers in different depositions -- exactly what
    the pre-SIFTS-fix ERalpha bug looked like, and exactly what this
    display exists to surface instead of silently hiding."""
    index = {500: {'1AAA': 'A353', '1BBB': 'A268'}}
    result = format_native_numbering(index, 500)
    assert 'A268' in result and 'A353' in result
    assert 'INCONSISTENT' in result


def test_format_native_numbering_unknown_position_returns_placeholder():
    assert format_native_numbering({}, 999) == '?'


def test_residue_numbering_table_augments_without_mutating_input():
    mined_rows = [
        {'residue': 'ASP194.A', 'interaction': 'HBDonor', 'frequency': 0.8},
        {'residue': 'GLU353.A', 'interaction': 'HBAcceptor', 'frequency': 0.5},
    ]
    structures = [
        _structure('1AAA', {('A', 194): 194}),
        _structure('1BBB', {('A', 353): 353}),
    ]
    original_first_row = dict(mined_rows[0])

    out = residue_numbering_table(mined_rows, structures, top_n=None)

    assert out[0]['native_numbering'] == 'A194'
    assert out[1]['native_numbering'] == 'A353'
    assert mined_rows[0] == original_first_row  # input untouched


def test_residue_numbering_table_respects_top_n():
    mined_rows = [{'residue': f'GLY{i}.A', 'interaction': 'Hydrophobic', 'frequency': 1.0 - i * 0.01}
                  for i in range(20)]
    out = residue_numbering_table(mined_rows, structures=[], top_n=5)
    assert len(out) == 5
    assert all(r['native_numbering'] == '?' for r in out)  # no residue_map data supplied
