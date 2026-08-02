import numpy as np
from rdkit import Chem
from posegate.receptor_prep import load_receptor_mol

lig = Chem.MolFromMolFile('data/five_family/3IEO_ligand_h.sdf', removeHs=False)
rec = load_receptor_mol('data/five_family/3IEO_receptor_h.pkl')
lig_coords = lig.GetConformer().GetPositions()
rec_coords = rec.GetConformer().GetPositions()

ne2_idx, h1382, h1383 = 1381, 1382, 1383
o7_idx = 7

pos_d = rec_coords[ne2_idx]
pos_o = lig_coords[o7_idx]
dist_da = np.linalg.norm(pos_d - pos_o)
print('D(NE2)...A(O7) heavy-atom distance:', dist_da)

for h_idx in (h1382, h1383):
    pos_h = rec_coords[h_idx]
    dist_ha = np.linalg.norm(pos_h - pos_o)
    vec_hd = pos_d - pos_h
    vec_ha = pos_o - pos_h
    cos_angle = np.dot(vec_hd, vec_ha) / (np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha))
    angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
    print(f'H{h_idx}: H...A dist={dist_ha:.2f}  D-H-A angle={angle:.1f}')
