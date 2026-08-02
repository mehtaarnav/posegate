from rdkit import Chem

lig = Chem.MolFromMolFile('data/five_family/3IEO_ligand_h.sdf', removeHs=False)
donor_smarts = Chem.MolFromSmarts('[$([O,S,#7;+0]),$([Nv4+1]),$([n+]c[nH])]-[H]')
print('SMARTS parsed:', donor_smarts is not None)
matches = lig.GetSubstructMatches(donor_smarts)
print('donor matches:', matches)
for idx in set(m[0] for m in matches):
    print('  donor atom', idx, lig.GetAtomWithIdx(idx).GetSymbol())

# Also check acceptor pattern near residue 92 side - print all N/O with H neighbor explicitly
for atom in lig.GetAtoms():
    if atom.GetSymbol() in ('N', 'O'):
        h_neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetSymbol() == 'H']
        print(atom.GetIdx(), atom.GetSymbol(), 'H neighbors:', h_neighbors, 'formal charge:', atom.GetFormalCharge())
