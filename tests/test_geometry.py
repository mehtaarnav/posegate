# posegate/tests/test_geometry.py
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from posegate.autopsy import (
    get_vdw_radius,
    find_hydrogen_bonds,
    find_aromatic_contacts,
    evaluate_steric_clashes
)

def test_vdw_radius_lookup():
    assert get_vdw_radius('C') == 1.70
    assert get_vdw_radius('O') == 1.52
    assert get_vdw_radius('Unknown') == 1.70  # Default fallback

def test_hydrogen_bonds():
    """Methanol (donor) placed to H-bond with water (acceptor): water's O sits
    ~2.0 A past methanol's O-H bond, roughly co-linear, so ProLIF's HBDonor
    geometry criteria (distance + D-H...A angle) should reliably fire."""
    methanol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    AllChem.EmbedMolecule(methanol, randomSeed=1)

    water = Chem.AddHs(Chem.MolFromSmiles('O'))
    AllChem.EmbedMolecule(water, randomSeed=1)

    donor_o_idx = [a.GetIdx() for a in methanol.GetAtoms() if a.GetSymbol() == 'O'][0]
    # The H bonded to the oxygen specifically (methanol also has 3 methyl
    # H's, which all have degree==1 too, so that alone can't distinguish them).
    donor_h_idx = [n.GetIdx() for n in methanol.GetAtomWithIdx(donor_o_idx).GetNeighbors() if n.GetSymbol() == 'H'][0]
    acceptor_o_idx = [a.GetIdx() for a in water.GetAtoms() if a.GetSymbol() == 'O'][0]

    conf = methanol.GetConformer()
    pos_o = np.array(conf.GetAtomPosition(donor_o_idx))
    pos_h = np.array(conf.GetAtomPosition(donor_h_idx))

    # Place water's O ~2.0 A beyond the H, along the O-H bond direction.
    direction = (pos_h - pos_o) / np.linalg.norm(pos_h - pos_o)
    target_pos_a = pos_h + direction * 2.0

    water_conf = water.GetConformer()
    water_o_pos = np.array(water_conf.GetAtomPosition(acceptor_o_idx))
    shift = target_pos_a - water_o_pos
    for atom in water.GetAtoms():
        p = water_conf.GetAtomPosition(atom.GetIdx())
        water_conf.SetAtomPosition(atom.GetIdx(), (p.x + shift[0], p.y + shift[1], p.z + shift[2]))

    hbonds = find_hydrogen_bonds(methanol, water)

    assert len(hbonds) > 0
    assert any(h['type'] == 'L_Donor -> R_Acceptor' for h in hbonds)
    hit = next(h for h in hbonds if h['type'] == 'L_Donor -> R_Acceptor')
    assert 0.0 < hit['distance_A'] < 3.5
    assert hit['angle_deg'] > 120.0

def test_aromatic_contacts():
    """Two benzene rings stacked face-to-face ~3.8 A apart should register
    a FaceToFace pi-stacking contact."""
    ring_coords = []
    radius = 1.39  # aromatic C-C bond length puts atoms on this ring radius
    for i in range(6):
        angle = i * (np.pi / 3)
        ring_coords.append((radius * np.cos(angle), radius * np.sin(angle), 0.0))

    def make_benzene(z_offset):
        mol = Chem.MolFromSmiles('c1ccccc1')
        mol = Chem.RWMol(mol)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(ring_coords):
            conf.SetAtomPosition(i, (x, y, z + z_offset))
        mol.AddConformer(conf)
        return mol.GetMol()

    lig = make_benzene(z_offset=0.0)
    rec = make_benzene(z_offset=3.8)

    aromatic = find_aromatic_contacts(lig, rec)

    assert len(aromatic) > 0
    assert any(a['geometry'] in ('FaceToFace', 'EdgeToFace') for a in aromatic)

def test_steric_clashes():
    """Test that overlapping atoms are correctly flagged."""
    mol1 = Chem.AddHs(Chem.MolFromSmiles('C'))
    AllChem.EmbedMolecule(mol1, randomSeed=1)

    mol2 = Chem.AddHs(Chem.MolFromSmiles('C'))
    AllChem.EmbedMolecule(mol2, randomSeed=1)

    # Force a severe clash by placing mol2's carbon exactly on top of mol1's carbon
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    c1_pos = conf1.GetAtomPosition(0)
    conf2.SetAtomPosition(0, c1_pos)  # Overlap!

    clashes = evaluate_steric_clashes(mol1, mol2, threshold_overlap=0.1)

    # The two carbons should have an overlap of ~3.4 A (1.7 + 1.7 - 0.0)
    assert len(clashes) > 0
    assert clashes[0]['ligand_atom'] == 'C0'
    assert clashes[0]['receptor_atom'] == 'C0'
    assert clashes[0]['overlap_A'] >= 3.0
