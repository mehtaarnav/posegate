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


# Catalytic metal ions are part of the binding site, not solvent. PDBFixer's
# removeHeterogens deletes every non-standard residue including these, which
# for a metalloenzyme removes the thing the ligand actually binds: carbonic
# anhydrase inhibitors coordinate the active-site Zn at about 2 A, so docking
# into a zinc-stripped receptor is not an approximation but a different
# problem. Alkali ions are excluded: Na and K in a crystal are almost always
# buffer, not catalytic.
METAL_IONS = {'ZN', 'MG', 'MN', 'FE', 'FE2', 'CU', 'CU1', 'CO', 'NI', 'CD', 'CA'}


def find_discarded_heterogens(pdb_path: str) -> set:
    """Names of HETATM residues that PDBFixer's removeHeterogens will
    delete and that this module does not put back (everything except
    water and the retained metal ions).

    A ligand that binds by coordinating a covalently- or non-covalently-
    bound cofactor -- heme iron, an FAD/NAD cofactor, ATP in a kinase's
    second site -- would be docked into a receptor missing that cofactor,
    with nothing to say so. This cannot be fixed generically the way
    metal retention was (a cofactor needs bonds and possibly missing-atom
    completion, not just re-inserted coordinates), so it is surfaced as a
    warning rather than silently handled.
    """
    found = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            resname = line[17:20].strip()
            if resname != 'HOH' and resname not in METAL_IONS:
                found.add(resname)
    return found


def warn_on_discarded_heterogens(pdb_path: str) -> None:
    discarded = find_discarded_heterogens(pdb_path)
    if discarded:
        print(f"receptor_prep: WARNING - removing non-metal heterogen(s) {sorted(discarded)} "
              f"from {pdb_path}. If any of these is a cofactor that forms part of the binding "
              f"site (e.g. a heme, ATP, NAD or FAD), the prepared receptor will be missing it "
              f"and docked poses will not reflect its presence.")


def pdb_atom_name(name: str) -> str:
    """Formats an atom name into PDB's 4-column field without truncating it.

    The previous formatting, f' {name:<3}'[:4], pads to 4 characters by
    prepending a space and then truncates to 4, which drops the last
    character of any 4-character name. HD11 and HD12 (the two delta
    hydrogens of LEU/VAL-family residues) both became "HD1", making them
    indistinguishable by name. Nothing currently reads AtomPDBResidueInfo's
    name field, so no existing result depends on this, but any output that
    round-trips through it (e.g. a written-out PDB) would silently lose the
    distinction.
    """
    return name if len(name) >= 4 else f' {name:<3}'


def extract_metal_ions(pdb_path: str):
    """Reads catalytic metal ions out of a PDB, before PDBFixer discards them.

    Returns a list of (resname, resnum, chain, element, (x, y, z)).
    """
    ions = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            resname = line[17:20].strip()
            if resname not in METAL_IONS:
                continue
            element = line[76:78].strip() or resname
            ions.append((
                resname,
                int(line[22:26]),
                line[21],
                element.capitalize(),
                (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            ))
    return ions


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
    # Read the metals out before PDBFixer removes them, and re-add them
    # below as isolated atoms. They carry no bonds: the coordination
    # geometry is what matters to interaction detection, and inventing
    # covalent bonds to a metal would be worse than leaving it unbonded.
    warn_on_discarded_heterogens(pdb_path)
    metal_ions = extract_metal_ions(pdb_path)

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
        info.SetName(pdb_atom_name(atom.name))
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

    metal_indices = []
    for resname, resnum, chain, element, _ in metal_ions:
        a = Chem.Atom(element)
        a.SetNoImplicit(True)  # a bare ion must not acquire implicit hydrogens
        info = Chem.AtomPDBResidueInfo()
        info.SetName(pdb_atom_name(resname))
        info.SetResidueName(resname)
        info.SetResidueNumber(resnum)
        info.SetChainId(chain)
        a.SetMonomerInfo(info)
        metal_indices.append(mol.AddAtom(a))
    if metal_ions:
        listed = ', '.join(f"{r}{n}.{c}" for r, n, c, _, _ in metal_ions)
        print(f"receptor_prep: kept {len(metal_ions)} metal ion(s) in the receptor: {listed}")

    conf = Chem.Conformer(mol.GetNumAtoms())
    positions = fixer.positions.value_in_unit(fixer.positions.unit)
    for omm_idx, rd_idx in atom_map.items():
        x, y, z = positions[omm_idx]
        conf.SetAtomPosition(rd_idx, Point3D(x * 10, y * 10, z * 10))  # nm -> Angstrom
    for rd_idx, (_, _, _, _, (x, y, z)) in zip(metal_indices, metal_ions):
        conf.SetAtomPosition(rd_idx, Point3D(x, y, z))  # already Angstrom
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
    warn_on_discarded_heterogens(pdb_path)
    metal_ions = extract_metal_ions(pdb_path)

    fixer = PDBFixer(filename=pdb_path)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(out_path, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    # Re-insert the metals PDBFixer stripped, so the PDBQT handed to Vina
    # still contains them. Without this the ligand would be docked into an
    # apo site for any metalloenzyme.
    if metal_ions:
        with open(out_path) as f:
            lines = [l for l in f if not l.startswith(('END', 'CONECT'))]
        serial = 90000
        for resname, resnum, chain, element, (x, y, z) in metal_ions:
            serial += 1
            lines.append(
                f"HETATM{serial:>5} {resname:>4}{'':1}{resname:>3} {chain}{resnum:>4}{'':4}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}{'':10}{element.upper():>2}\n"
            )
        lines.append('END\n')
        with open(out_path, 'w') as f:
            f.writelines(lines)
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
