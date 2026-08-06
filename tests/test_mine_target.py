# posegate/tests/test_mine_target.py
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from mine_target import detect_ligand_resname


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
