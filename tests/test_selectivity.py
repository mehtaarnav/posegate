# posegate/tests/test_selectivity.py
import json

import pytest

from posegate.selectivity import (
    build_alignment_map, top_residue_numbers, compare_isoforms, fetch_uniprot_sequence,
)


def _write_mined(tmp_path, name, rows):
    path = tmp_path / f"{name}_mined.json"
    path.write_text(json.dumps({'mined': rows}))
    return str(path)


def _row(residue, freq, interaction='HBDonor'):
    return {'residue': residue, 'interaction': interaction, 'frequency': freq,
            'n_structures': 1, 'n_ensemble': 1, 'ci95': [0.0, 1.0]}


def test_top_residue_numbers_ranks_by_frequency_and_dedupes(tmp_path):
    path = _write_mined(tmp_path, "s", [
        _row('GLU106.A', 0.9),
        _row('GLU106.A', 0.9, interaction='Hydrophobic'),  # same residue, second interaction row
        _row('ASP83.A', 0.5),
        _row('LEU20.A', 0.3, interaction='VdWContact'),  # excluded
    ])
    assert top_residue_numbers(path, top_n=10) == [106, 83]


def test_top_residue_numbers_respects_top_n(tmp_path):
    path = _write_mined(tmp_path, "s", [_row(f'GLY{i}.A', 1.0 - i * 0.01) for i in range(20)])
    assert len(top_residue_numbers(path, top_n=5)) == 5


def test_top_residue_numbers_missing_file_raises_actionable_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="mine_target.py"):
        top_residue_numbers(str(tmp_path / "nope.json"), top_n=10)


def test_top_residue_numbers_empty_mined_list_returns_empty(tmp_path):
    path = _write_mined(tmp_path, "s", [])
    assert top_residue_numbers(path, top_n=10) == []


def test_top_residue_numbers_rejects_bad_top_n(tmp_path):
    path = _write_mined(tmp_path, "s", [_row('GLU106.A', 0.9)])
    with pytest.raises(ValueError):
        top_residue_numbers(path, top_n=0)


def test_build_alignment_map_identity_alignment_maps_every_position():
    seq = "MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLRILNNGHAFNVEFDDSQDKAVLKGGPLDGTYRLIQFHFHWGSLDGQGSEHTVDKKKYAAELHLVHWNTKYGDFGKAVQQPDGLAVLGIFLKVGSAKPGLQKVVDVLDSIKTKGKSADFTNFDPRGLLPESLDYWTYPGSLTTPPLLECVTWIVLKEPISVSSEQVLKFRKLNFNGEGEPEELMVDNWRPAQPLKNRQIKASFK"
    m = build_alignment_map(seq, seq)
    assert m[1] == 1
    assert m[len(seq)] == len(seq)
    assert len(m) == len(seq)


def test_build_alignment_map_finds_correct_offset_shift():
    """A sequence embedded with a prefix insertion (simulating an
    isoform with an extra N-terminal domain, like CA9's PG domain) must
    still map corresponding residues to the same relative position, not
    the same absolute number."""
    core = "MSHHWGYGKHNGPEHWHKDFPIAKGERQSPVDIDTHTAKYDPSLKPLSVSYDQATSLR"
    shifted = "XXXXXXXXXX" + core  # 10-residue N-terminal insertion
    m = build_alignment_map(core, shifted)
    # position 1 in `core` should map to position 11 in `shifted` (after the insertion)
    assert m[1] == 11
    assert m[len(core)] == len(shifted)


def test_build_alignment_map_rejects_empty_sequence():
    with pytest.raises(ValueError):
        build_alignment_map("", "ACDEFG")


def test_compare_isoforms_requires_reference_present(monkeypatch, tmp_path):
    isoforms = {'A': {'acc': 'X', 'mined': str(tmp_path / "a.json")}}
    with pytest.raises(ValueError, match="not in isoforms"):
        compare_isoforms(isoforms, reference='B')


def test_compare_isoforms_requires_at_least_two(tmp_path):
    isoforms = {'A': {'acc': 'X', 'mined': str(tmp_path / "a.json")}}
    with pytest.raises(ValueError, match="at least 2"):
        compare_isoforms(isoforms, reference='A')


def test_compare_isoforms_shared_pairwise_and_unique_categories(monkeypatch, tmp_path):
    """Three synthetic isoforms sharing one core sequence exactly (so
    alignment is trivial identity), with engineered top-N residue sets
    that exercise every category compare_isoforms reports: shared by
    all three, shared by exactly one pair, and unique to one isoform."""
    seq = "ACDEFGHIKLMNPQRSTVWY" * 3  # 60 residues, identical across all three
    monkeypatch.setattr('posegate.selectivity.fetch_uniprot_sequence', lambda acc: seq)

    mined = {
        'A': _write_mined(tmp_path, "A", [_row('X5.A', 1.0), _row('X10.A', 0.9), _row('X15.A', 0.5)]),
        'B': _write_mined(tmp_path, "B", [_row('X5.A', 1.0), _row('X10.A', 0.9), _row('X20.A', 0.4)]),
        'C': _write_mined(tmp_path, "C", [_row('X5.A', 1.0), _row('X30.A', 0.6)]),
    }
    isoforms = {n: {'acc': f'acc{n}', 'mined': mined[n]} for n in ('A', 'B', 'C')}

    result = compare_isoforms(isoforms, reference='A', top_n=10)

    assert result['shared_by_all'] == {5}
    assert result['pairwise_only'][('A', 'B')] == {10}
    assert result['unique_to']['B'] == {20}
    assert result['unique_to']['C'] == {30}
    assert result['unique_to']['A'] == {15}


def test_fetch_uniprot_sequence_raises_clear_error_on_bad_response(monkeypatch):
    class _FakeResp:
        text = "not fasta at all"
        def raise_for_status(self):
            pass

    monkeypatch.setattr('requests.get', lambda *a, **k: _FakeResp())
    with pytest.raises(ValueError, match="not FASTA"):
        fetch_uniprot_sequence("NOPE")
