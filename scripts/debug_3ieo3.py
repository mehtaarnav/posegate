from rdkit import Chem

lig = Chem.MolFromMolFile('data/five_family/3IEO_ligand_h.sdf', removeHs=False)
acceptor_smarts = Chem.MolFromSmarts(
    '[$([N&!$([NX3]-*=[O,N,P,S])&!$([ND2v3^2+0](-[H])-[CD4v4H1^3]-[CD2^2+0]=O)&!$([NX3]-[a])&!$([Nv4+1])&!$(N=C(-[C,N])-N)]),'
    '$([n+0&!X3&!$([n&r5]:[n+&r5])]),'
    '$([O&!$([OX2](C)C=O)&!$(O(~a)~a)&!$(O=N-*)&!$([O-]-N=O)]),'
    '$([o+0]),'
    '$([F&$(F-[#6])&!$(F-[#6][F,Cl,Br,I])])]'
)
print('parsed:', acceptor_smarts is not None)
matches = set(m[0] for m in lig.GetSubstructMatches(acceptor_smarts))
print('acceptor atom indices:', matches)

atom7 = lig.GetAtomWithIdx(7)
print('atom 7:', atom7.GetSymbol(), 'neighbors:', [(n.GetIdx(), n.GetSymbol()) for n in atom7.GetNeighbors()])
print('atom 7 is acceptor match:', 7 in matches)

# print full connectivity around atom 7 to see the ester/ether context
for n in atom7.GetNeighbors():
    print('  neighbor', n.GetIdx(), n.GetSymbol(), 'its neighbors:', [(nn.GetIdx(), nn.GetSymbol(), lig.GetBondBetweenAtoms(n.GetIdx(), nn.GetIdx()).GetBondTypeAsDouble()) for nn in n.GetNeighbors()])
