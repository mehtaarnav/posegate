# posegate/tests/test_conserved_contacts.py
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from posegate.conserved_contacts import mine_conserved_contacts

def _make_structure(tmp_path, name, shift):
    """A tiny synthetic 'structure': a methanol ligand H-bonding to a water
    receptor placed at a fixed offset, so the same (residue, interaction)
    pair should be mined as conserved across repeated calls with the same
    shift, self-contained (no downloaded PDB data required)."""
    methanol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    AllChem.EmbedMolecule(methanol, randomSeed=1)
    donor_o_idx = 1
    donor_h_idx = [n.GetIdx() for n in methanol.GetAtomWithIdx(donor_o_idx).GetNeighbors() if n.GetSymbol() == 'H'][0]

    water = Chem.AddHs(Chem.MolFromSmiles('O'))
    AllChem.EmbedMolecule(water, randomSeed=1)

    conf = methanol.GetConformer()
    import numpy as np
    pos_o = np.array(conf.GetAtomPosition(donor_o_idx))
    pos_h = np.array(conf.GetAtomPosition(donor_h_idx))
    direction = (pos_h - pos_o) / np.linalg.norm(pos_h - pos_o)
    target_pos_a = pos_h + direction * 2.0

    water_conf = water.GetConformer()
    water_o_pos = np.array(water_conf.GetAtomPosition(0))
    move = target_pos_a - water_o_pos + np.array(shift)
    for atom in water.GetAtoms():
        p = water_conf.GetAtomPosition(atom.GetIdx())
        water_conf.SetAtomPosition(atom.GetIdx(), (p.x + move[0], p.y + move[1], p.z + move[2]))

    lig_path = tmp_path / f"{name}_lig.sdf"
    rec_path = tmp_path / f"{name}_rec.pdb"
    Chem.MolToMolFile(methanol, str(lig_path))
    Chem.MolToPDBFile(water, str(rec_path))
    return {'pdb_id': name, 'ligand_sdf': str(lig_path), 'receptor_pdb': str(rec_path)}

def test_mine_conserved_contacts_finds_repeated_interaction(tmp_path):
    structures = [_make_structure(tmp_path, f"s{i}", shift=(0, 0, 0)) for i in range(3)]
    results = mine_conserved_contacts(structures)

    assert len(results) > 0
    top = results[0]
    assert top['frequency'] == 1.0
    assert top['n_structures'] == 3

def test_mine_conserved_contacts_reports_partial_frequency(tmp_path):
    # Two structures with the same water position (should H-bond), one with
    # the water pushed far away (should not) -> that interaction should
    # show up at 2/3 frequency, not 3/3.
    structures = [
        _make_structure(tmp_path, "near1", shift=(0, 0, 0)),
        _make_structure(tmp_path, "near2", shift=(0, 0, 0)),
        _make_structure(tmp_path, "far", shift=(20, 20, 20)),
    ]
    results = mine_conserved_contacts(structures)

    hbond_hits = [r for r in results if 'HB' in r['interaction']]
    assert len(hbond_hits) > 0
    assert hbond_hits[0]['n_structures'] == 2
    assert hbond_hits[0]['frequency'] == pytest.approx(2 / 3, abs=0.01)
