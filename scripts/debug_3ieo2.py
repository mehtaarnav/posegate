import numpy as np
from rdkit import Chem
from posegate.receptor_prep import load_receptor_mol

lig = Chem.MolFromMolFile('data/five_family/3IEO_ligand_h.sdf', removeHs=False)
rec = load_receptor_mol('data/five_family/3IEO_receptor_h.pkl')

lig_coords = lig.GetConformer().GetPositions()
rec_coords = rec.GetConformer().GetPositions()

# Gln92 sidechain amide: NE2-HE21/HE22 (donor), OE1 (acceptor)
gln92_atoms = [(a.GetIdx(), a.GetPDBResidueInfo().GetName().strip(), a.GetSymbol())
               for a in rec.GetAtoms() if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetResidueNumber() == 92]
print('Gln92 atoms:', gln92_atoms)

ne2_idx = [i for i, n, s in gln92_atoms if n == 'NE2'][0]
h_idxs = [n.GetIdx() for n in rec.GetAtomWithIdx(ne2_idx).GetNeighbors() if n.GetSymbol() == 'H']
print('NE2 H neighbors:', h_idxs)

# ligand O atoms
lig_o_idxs = [a.GetIdx() for a in lig.GetAtoms() if a.GetSymbol() in ('O', 'N')]
for h_idx in h_idxs:
    h_pos = rec_coords[h_idx]
    for o_idx in lig_o_idxs:
        o_pos = lig_coords[o_idx]
        dist = np.linalg.norm(h_pos - o_pos)
        if dist < 4.0:
            print(f'  H{h_idx} - ligand atom {o_idx} ({lig.GetAtomWithIdx(o_idx).GetSymbol()}): {dist:.2f} A')
