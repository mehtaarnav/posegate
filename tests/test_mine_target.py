# posegate/tests/test_mine_target.py
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
import mine_target
from mine_target import detect_ligand_resname, global_pdb_prevalence, PREVALENCE_THRESHOLD


def _hetatm_line(resname, chain, resnum, n_atoms_marker=1):
    # Column-exact PDB HETATM format: resName at [17:20], chainID at [21],
    # resSeq at [22:26] -- must match what detect_ligand_resname reads.
    return (f"HETATM{n_atoms_marker:>5}  C   {resname:>3} {chain}{resnum:>4}      "
            f"0.000   0.000   0.000  1.00  0.00           C\n")


def test_detect_ligand_resname_prefers_real_ligand_over_larger_cofactor(tmp_path):
    """Heme (and other prosthetic-group cofactors) must never be picked
    as 'the ligand' even though they have more atoms than most drug-like
    ligands -- found on COX-1, a heme-dependent peroxidase, where HEM
    outsizing the real bound NSAID crashed the whole ensemble."""
    pdb = tmp_path / "cox_like.pdb"
    lines = []
    # HEM: 43 atoms, larger than the real ligand -- must still lose.
    for i in range(43):
        lines.append(_hetatm_line('HEM', 'A', 901, n_atoms_marker=i + 1))
    # Real ligand (e.g. ibuprofen-like): fewer atoms.
    for i in range(15):
        lines.append(_hetatm_line('IBP', 'A', 902, n_atoms_marker=100 + i))
    pdb.write_text(''.join(lines))

    assert detect_ligand_resname(str(pdb)) == 'IBP'


def test_detect_ligand_resname_prefers_real_ligand_over_glycosylation_sugar(tmp_path):
    """N-acetylglucosamine (NAG) and other glycosylation sugars are
    present on essentially every glycoprotein and must never be picked
    as 'the ligand' -- found on ERalpha and again on COX-2 (see
    conversation), where NAG outcompeted the real bound inhibitor."""
    pdb = tmp_path / "glyco_like.pdb"
    lines = []
    for i in range(20):
        lines.append(_hetatm_line('NAG', 'A', 601, n_atoms_marker=i + 1))
    for i in range(12):
        lines.append(_hetatm_line('SAL', 'A', 602, n_atoms_marker=100 + i))
    pdb.write_text(''.join(lines))

    assert detect_ligand_resname(str(pdb)) == 'SAL'


def test_detect_ligand_resname_returns_none_when_only_cofactors_present(tmp_path):
    pdb = tmp_path / "apo_like.pdb"
    lines = [_hetatm_line('HEM', 'A', 901, n_atoms_marker=i + 1) for i in range(43)]
    pdb.write_text(''.join(lines))

    assert detect_ligand_resname(str(pdb)) is None


def test_detect_ligand_resname_excludes_common_component_not_in_static_list(tmp_path, monkeypatch):
    """The whole point of the prevalence check: a contaminant class that
    was NEVER manually added to LIKELY_NON_LIGAND must still be excluded,
    because it's globally common across the PDB -- this is what makes
    the fix generalize instead of requiring a new hand-added entry for
    every future cofactor/buffer/detergent class, the way HEM, NAG and
    BOG each did before this existed (see conversation)."""
    mine_target._prevalence_cache.clear()
    pdb = tmp_path / "common_additive.pdb"
    lines = []
    # 'ZZZ' stands in for an unlisted-but-common additive: more atoms
    # than the real ligand, NOT in LIKELY_NON_LIGAND.
    for i in range(30):
        lines.append(_hetatm_line('ZZZ', 'A', 701, n_atoms_marker=i + 1))
    for i in range(12):
        lines.append(_hetatm_line('RGL', 'A', 702, n_atoms_marker=100 + i))
    pdb.write_text(''.join(lines))

    def fake_prevalence(comp_id, timeout=15):
        return {'ZZZ': 9000, 'RGL': 8}[comp_id]
    monkeypatch.setattr(mine_target, 'global_pdb_prevalence', fake_prevalence)

    assert detect_ligand_resname(str(pdb)) == 'RGL'


def test_global_pdb_prevalence_caches_results(monkeypatch):
    mine_target._prevalence_cache.clear()
    call_count = {'n': 0}

    class _FakeResp:
        status_code = 200
        def raise_for_status(self):
            pass
        def json(self):
            return {'total_count': 42}

    def fake_post(*a, **k):
        call_count['n'] += 1
        return _FakeResp()

    monkeypatch.setattr(mine_target.requests, 'post', fake_post)
    assert global_pdb_prevalence('XYZ') == 42
    assert global_pdb_prevalence('XYZ') == 42
    assert call_count['n'] == 1


def test_global_pdb_prevalence_fails_open_on_network_error(monkeypatch):
    import requests
    mine_target._prevalence_cache.clear()

    def fake_post(*a, **k):
        raise requests.RequestException("network down")
    monkeypatch.setattr(mine_target.requests, 'post', fake_post)

    assert global_pdb_prevalence('ABC') == -1


def test_detect_ligand_resname_accepts_candidate_when_prevalence_check_fails(tmp_path, monkeypatch):
    """A network hiccup during the prevalence check must not silently
    drop a real structure from the ensemble -- fail open, not closed."""
    mine_target._prevalence_cache.clear()
    pdb = tmp_path / "network_down.pdb"
    pdb.write_text(''.join(_hetatm_line('QRS', 'A', 501, n_atoms_marker=i + 1) for i in range(10)))

    monkeypatch.setattr(mine_target, 'global_pdb_prevalence', lambda c, timeout=15: -1)

    assert detect_ligand_resname(str(pdb)) == 'QRS'
