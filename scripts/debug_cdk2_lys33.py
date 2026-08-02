import json
import numpy as np
from rdkit import Chem
from posegate.receptor_prep import load_receptor_mol

with open('data/ensemble_cdk2/cdk2_prepped.json') as f:
    structures = json.load(f)

for s in structures:
    lig = Chem.MolFromMolFile(s['ligand_sdf'], removeHs=False)
    rec = load_receptor_mol(s['receptor_pdb'])
    if lig is None or rec is None:
        print(s['pdb_id'], 'FAILED TO LOAD')
        continue

    lys33_atoms = [a for a in rec.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueNumber() == 33]
    nz = [a for a in lys33_atoms if a.GetPDBResidueInfo().GetName().strip() == 'NZ']
    if not nz:
        print(s['pdb_id'], 'no NZ atom found on Lys33')
        continue
    nz_idx = nz[0].GetIdx()
    h_idxs = [n.GetIdx() for n in rec.GetAtomWithIdx(nz_idx).GetNeighbors() if n.GetSymbol() == 'H']

    rec_coords = rec.GetConformer().GetPositions()
    lig_coords = lig.GetConformer().GetPositions()
    lig_hetero = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() in ('N', 'O')]

    best = None
    for h_idx in h_idxs:
        pos_h = rec_coords[h_idx]
        pos_d = rec_coords[nz_idx]
        for o_idx in lig_hetero:
            pos_a = lig_coords[o_idx]
            dist_ha = np.linalg.norm(pos_h - pos_a)
            if dist_ha > 3.0:
                continue
            vec_hd = pos_d - pos_h
            vec_ha = pos_a - pos_h
            cos_angle = np.dot(vec_hd, vec_ha) / (np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha))
            angle = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
            dist_da = np.linalg.norm(pos_d - pos_a)
            if best is None or dist_ha < best[0]:
                best = (dist_ha, dist_da, angle, o_idx)

    if best:
        print(f"{s['pdb_id']}: closest Lys33-NZ...ligand contact: H...A={best[0]:.2f} D...A={best[1]:.2f} D-H-A angle={best[2]:.1f}")
    else:
        print(f"{s['pdb_id']}: no Lys33-NZ contact within 3.0 A of any ligand N/O")
