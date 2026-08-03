# posegate/posegate/receptor_prep.py
"""Prepares a raw receptor PDB into a properly-bonded RDKit molecule.

Loading a receptor via Chem.MolFromPDBFile alone is unreliable for full
proteins: RDKit's PDB parser either falls back to naive proximity-based
bonding (which can create spurious bonds between spatially close but
unconnected atoms, e.g. residue i's carbonyl C to residue i+1's carbonyl
O on a tight turn) or, with that fallback disabled, its residue-template
matcher can silently fail to bond the vast majority of a multi-residue
chain (observed: only the first 1-2 residues bonded, the other ~297
completely unbonded, on an otherwise perfectly ordinary CDK2 structure).

PDBFixer/OpenMM's Topology object already computes correct, chemically
valid bonds during its own residue-template matching (that's how OpenMM
runs simulations) via `PDBFixer.topology.bonds()`. This module builds an
RDKit Mol directly from that Topology + positions instead of writing PDB
text and re-guessing bonds from it, which is both what throws the
already-correct bonds away and what every native RDKit/MDAnalysis PDB
bond-guesser choked on above.
"""

import pickle

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Geometry import Point3D


# Ring atoms of the aromatic side chains, by PDB atom name, in ring
# connectivity order. Tryptophan contributes two fused rings and appears
# twice; the shared CD2-CE2 bond is set aromatic by whichever comes first.
AROMATIC_SIDECHAIN_RINGS = {
    'PHE': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    'TYR': [['CG', 'CD1', 'CE1', 'CZ', 'CE2', 'CD2']],
    'HIS': [['CG', 'ND1', 'CE1', 'NE2', 'CD2']],
    'TRP': [['CG', 'CD1', 'NE1', 'CE2', 'CD2'],
            ['CD2', 'CE2', 'CZ2', 'CH2', 'CZ3', 'CE3']],
}


def assign_sidechain_aromaticity(mol: Chem.RWMol) -> int:
    """Flags the aromatic side-chain rings of PHE, TYR, HIS and TRP.

    Bonds taken from OpenMM's Topology carry no bond order, so every bond
    is created single and nothing in the molecule is aromatic. ProLIF
    detects pi-stacking (FaceToFace/EdgeToFace) from aromatic ring
    perception, so without this no aromatic interaction can ever be
    reported, on any structure: the feature silently reads zero rather
    than failing. Ring membership for the standard residues is known from
    their names, so it is assigned from a template rather than inferred.

    Returns the number of rings flagged.
    """
    by_residue = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        resname = info.GetResidueName().strip()
        if resname not in AROMATIC_SIDECHAIN_RINGS:
            continue
        key = (resname, info.GetResidueNumber(), info.GetChainId())
        by_residue.setdefault(key, {})[info.GetName().strip()] = atom.GetIdx()

    flagged = 0
    for (resname, _, _), atoms in by_residue.items():
        for ring in AROMATIC_SIDECHAIN_RINGS[resname]:
            if not all(name in atoms for name in ring):
                continue  # incomplete side chain; leave it alone
            idxs = [atoms[name] for name in ring]
            ok = True
            for i in range(len(idxs)):
                bond = mol.GetBondBetweenAtoms(idxs[i], idxs[(i + 1) % len(idxs)])
                if bond is None:
                    ok = False
                    break
            if not ok:
                continue
            for i in range(len(idxs)):
                mol.GetBondBetweenAtoms(
                    idxs[i], idxs[(i + 1) % len(idxs)]
                ).SetBondType(Chem.BondType.AROMATIC)
                mol.GetAtomWithIdx(idxs[i]).SetIsAromatic(True)
            flagged += 1
    return flagged


def prepare_receptor_mol(pdb_path: str) -> Chem.Mol:
    """Cleans heterogens/adds missing atoms+hydrogens via PDBFixer, then
    builds an RDKit Mol directly from PDBFixer's own Topology bonds."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    mol = Chem.RWMol()
    atom_map = {}
    skipped_residues = set()
    for atom in fixer.topology.atoms():
        # Residues numbered below 1 are cloning/expression-tag artifacts
        # (e.g. 6Q3F carries SER-3, PRO-2, GLU-1, PHE0 ahead of CDK2's
        # Met1). They are never part of a binding site, and ProLIF stores
        # residue numbers as uint32, so leaving them in raises
        # "OverflowError: Python integer -3 out of bounds for uint32" when
        # the fingerprint is built. Dropped here, and reported rather than
        # removed silently.
        if int(atom.residue.id) < 1:
            skipped_residues.add((atom.residue.name, int(atom.residue.id), atom.residue.chain.id))
            continue
        a = Chem.Atom(atom.element.symbol)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f' {atom.name:<3}'[:4])
        info.SetResidueName(atom.residue.name)
        info.SetResidueNumber(int(atom.residue.id))
        info.SetChainId(atom.residue.chain.id)
        a.SetMonomerInfo(info)
        idx = mol.AddAtom(a)
        atom_map[atom.index] = idx

    if skipped_residues:
        listed = ', '.join(f"{n}{i}.{c}" for n, i, c in sorted(skipped_residues, key=lambda r: r[1]))
        print(f"receptor_prep: dropped {len(skipped_residues)} residue(s) numbered below 1 "
              f"(expression-tag artifacts, not part of the binding site): {listed}")

    for bond in fixer.topology.bonds():
        # Bonds to a dropped tag residue have no counterpart here.
        if bond.atom1.index not in atom_map or bond.atom2.index not in atom_map:
            continue
        i, j = atom_map[bond.atom1.index], atom_map[bond.atom2.index]
        if not mol.GetBondBetweenAtoms(i, j):
            mol.AddBond(i, j, Chem.BondType.SINGLE)

    conf = Chem.Conformer(mol.GetNumAtoms())
    positions = fixer.positions.value_in_unit(fixer.positions.unit)
    for omm_idx, rd_idx in atom_map.items():
        x, y, z = positions[omm_idx]
        conf.SetAtomPosition(rd_idx, Point3D(x * 10, y * 10, z * 10))  # nm -> Angstrom
    mol.AddConformer(conf)

    # Topology.bonds() carries no bond order, so every bond above is
    # single and nothing is aromatic. Flag the aromatic side chains from
    # their residue templates before sanitizing: ProLIF derives
    # pi-stacking from aromatic ring perception, so skipping this makes
    # FaceToFace/EdgeToFace unreportable on every structure, which reads
    # as "no aromatic interactions found" rather than as an error.
    assign_sidechain_aromaticity(mol)

    final = mol.GetMol()
    # Kekulization and RDKit's own aromaticity perception stay off: bond
    # orders are unknown, so there is no Kekule structure to find, and
    # re-perceiving would discard the flags just assigned. Ring info is
    # still computed, which is what ProLIF needs alongside those flags.
    ops = Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_PROPERTIES
    Chem.SanitizeMol(final, sanitizeOps=ops)
    return final


def write_clean_receptor_pdb(pdb_path: str, out_path: str) -> str:
    """Writes a heterogen-free, hydrogenated receptor PDB.

    Docking must not see the co-crystallized ligand: the raw crystal file
    still contains it, and docking into that structure silently blocks the
    real binding site with a copy of the native ligand. Producing this file
    used to be a manual step, done once for BRD4 and never scripted, which
    left every new target unable to run batch_dock.py at all.

    Unlike prepare_receptor_pickle, the output here is PDB text and so does
    lose the computed bonds. That is fine for its two consumers: OpenBabel
    only needs atoms and elements to write a PDBQT for Vina, and PDBFixer
    recomputes its own Topology when this file is read back for the
    pickle. Anything needing correct bonds should use the pickle.
    """
    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(out_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    return out_path


def prepare_receptor_pickle(pdb_path: str, out_path: str) -> str:
    """Prepares a receptor and pickles the resulting RDKit Mol. Pickling
    (rather than writing PDB/SDF text) is the only lossless round trip:
    plain PDB text drops the bonds we just computed (standard residues
    aren't given explicit CONECT records by convention) and plain SDF/MOL
    text has no field for AtomPDBResidueInfo (residue name/number/chain),
    which posegate's per-residue conserved-contact analysis depends on."""
    mol = prepare_receptor_mol(pdb_path)
    with open(out_path, 'wb') as f:
        pickle.dump(mol, f)
    return out_path


def load_receptor_mol(pkl_path: str) -> Chem.Mol:
    with open(pkl_path, 'rb') as f:
        return pickle.load(f)
