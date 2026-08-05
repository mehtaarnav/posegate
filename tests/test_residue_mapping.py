# posegate/tests/test_residue_mapping.py
from posegate.residue_mapping import build_residue_map, remap_residue_numbers


class _FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


def _fake_sifts_payload(pdb_id, chain, label_start, label_end, unp_start):
    return {
        pdb_id: {
            'UniProt': {
                'P03372': {
                    'identifier': 'ESR1_HUMAN',
                    'mappings': [{
                        'chain_id': chain,
                        'start': {'author_residue_number': None, 'residue_number': label_start},
                        'end': {'author_residue_number': None, 'residue_number': label_end},
                        'unp_start': unp_start,
                        'unp_end': unp_start + (label_end - label_start),
                    }],
                }
            }
        }
    }


def _fake_residue_listing(pdb_id, chain, label_to_author):
    return {
        pdb_id: {'molecules': [{'chains': [{
            'struct_asym_id': chain,
            'residues': [{'residue_number': label, 'author_residue_number': auth,
                          'residue_name': 'ALA', 'author_insertion_code': ''}
                         for label, auth in label_to_author.items()],
        }]}]}
    }


def test_build_residue_map_applies_a_constant_offset_within_a_segment(monkeypatch):
    """Two structures whose author numbering for the same physical
    residues differs by a fixed offset (the ordinary case: different
    construct start points, no internal indel) must map to the *same*
    UniProt residue numbers, since that's the whole point -- making
    (chain, resnum) mean the same physical position across structures
    that were deposited with different author numbering."""
    import posegate.residue_mapping as rm

    def fake_get_a(url, timeout):
        if 'residue_listing' in url:
            return _FakeResponse(_fake_residue_listing('1aaa', 'A', {i: 300 + (i - 1) for i in range(1, 12)}))
        return _FakeResponse(_fake_sifts_payload('1aaa', 'A', label_start=1, label_end=11, unp_start=300))

    def fake_get_b(url, timeout):
        if 'residue_listing' in url:
            return _FakeResponse(_fake_residue_listing('1bbb', 'A', {i: i for i in range(1, 12)}))
        return _FakeResponse(_fake_sifts_payload('1bbb', 'A', label_start=1, label_end=11, unp_start=300))

    monkeypatch.setattr(rm.requests, 'get', fake_get_a)
    map_a = build_residue_map('1aaa', 'P03372')
    monkeypatch.setattr(rm.requests, 'get', fake_get_b)
    map_b = build_residue_map('1bbb', 'P03372')

    # Same physical residue (5th position into the segment) in each
    # structure's own author numbering -> identical UniProt number.
    assert map_a[('A', 305)] == map_b[('A', 6)] == 305


def test_build_residue_map_raises_when_accession_not_in_sifts(monkeypatch):
    import posegate.residue_mapping as rm
    monkeypatch.setattr(rm.requests, 'get', lambda url, timeout: _FakeResponse({'1ccc': {'UniProt': {}}}))
    try:
        build_residue_map('1ccc', 'P03372')
        assert False, "expected ValueError"
    except ValueError as e:
        assert 'P03372' in str(e)


def test_remap_residue_numbers_rewrites_only_mapped_atom_lines(tmp_path):
    pdb_in = tmp_path / "in.pdb"
    pdb_out = tmp_path / "out.pdb"
    pdb_in.write_text(
        "ATOM      1  N   GLU A 305      23.907  -6.371  21.286  1.00 24.42           N\n"
        "ATOM      2  CA  GLU A 305      23.055  -5.177  21.259  1.00 26.41           C\n"
        "HETATM    3  O   HOH B  10      10.000   1.000   2.000  1.00 30.00           O\n"
    )
    residue_map = {('A', 305): 353}

    n = remap_residue_numbers(str(pdb_in), str(pdb_out), residue_map)

    assert n == 2
    lines = pdb_out.read_text().splitlines()
    assert lines[0][22:26].strip() == '353'
    assert lines[1][22:26].strip() == '353'
    # unmapped chain/residue (the water) passes through unchanged
    assert lines[2][22:26].strip() == '10'
