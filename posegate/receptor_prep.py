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
    for atom in fixer.topology.atoms():
        a = Chem.Atom(atom.element.symbol)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f' {atom.name:<3}'[:4])
        info.SetResidueName(atom.residue.name)
        info.SetResidueNumber(int(atom.residue.id))
        info.SetChainId(atom.residue.chain.id)
        a.SetMonomerInfo(info)
        idx = mol.AddAtom(a)
        atom_map[atom.index] = idx

    for bond in fixer.topology.bonds():
        i, j = atom_map[bond.atom1.index], atom_map[bond.atom2.index]
        if not mol.GetBondBetweenAtoms(i, j):
            mol.AddBond(i, j, Chem.BondType.SINGLE)

    conf = Chem.Conformer(mol.GetNumAtoms())
    positions = fixer.positions.value_in_unit(fixer.positions.unit)
    for omm_idx, rd_idx in atom_map.items():
        x, y, z = positions[omm_idx]
        conf.SetAtomPosition(rd_idx, Point3D(x * 10, y * 10, z * 10))  # nm -> Angstrom
    mol.AddConformer(conf)

    final = mol.GetMol()
    # All bonds are added as single-order above (Topology.bonds() doesn't
    # reliably expose order for every residue), so skip kekulization and
    # aromaticity perception, which would otherwise reject this on
    # aromatic side chains; this doesn't affect the geometry/connectivity
    # ProLIF's interaction detection actually needs.
    ops = Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_SETAROMATICITY ^ Chem.SANITIZE_PROPERTIES
    Chem.SanitizeMol(final, sanitizeOps=ops)
    return final


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
