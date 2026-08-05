# posegate/posegate/autopsy.py
import prolif
from rdkit import Chem
from typing import List, Dict, Any, Optional, Sequence, Tuple, Union

from posegate.receptor_prep import load_receptor_mol

# A conserved-contact residue, as (3-letter name, number, chain id).
ConservedResidue = Tuple[str, int, str]

# BRD4's Asn140 acetyl-lysine-mimetic contact, the constraint this module
# was originally written around. Retained as the default only so existing
# BRD4 calls keep working; it is not a sensible default for any other
# target, and callers are expected to pass their own (from the miner).
DEFAULT_CONSERVED_RESIDUES: List[ConservedResidue] = [('ASN', 140, 'A')]

VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
    'S': 1.80, 'P': 1.80, 'F': 1.47, 'Cl': 1.75,
    'Br': 1.85, 'I': 1.98
}

IFP = Dict[Any, Dict[str, tuple]]

def get_vdw_radius(element: str) -> float:
    return VDW_RADII.get(element.upper(), 1.70)

def load_pose_mol(sdf_path: str) -> Optional[Chem.Mol]:
    """Loads a docked pose, tolerating output RDKit rejects outright.

    Poses reach us as SDF written by OpenBabel from Vina's PDBQT, and that
    round trip can produce valence or kekulization states RDKit's default
    parser refuses, returning None. Large peptidomimetic ligands hit this
    disproportionately, so treating an unparseable pose as a failed
    compound quietly biases a benchmark toward its smaller ligands.
    Parsing without sanitization and then sanitizing everything except the
    steps that fail on this output recovers the pose with its geometry,
    which is what interaction detection needs.
    """
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is not None:
        return mol
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE ^ Chem.SANITIZE_PROPERTIES,
        )
    except Exception:
        return None
    return mol


def build_ifp(ligand_mol: Chem.Mol, receptor_mol: Chem.Mol) -> IFP:
    """Computes the full ProLIF interaction fingerprint once per
    ligand/receptor pair; all interaction-detection functions below derive
    their result from this shared IFP rather than re-scanning atom pairs
    with hand-rolled geometry."""
    plig = prolif.Molecule.from_rdkit(ligand_mol, resname='LIG')
    prot = prolif.Molecule.from_rdkit(receptor_mol)
    fp = prolif.Fingerprint('all')
    return fp.generate(plig, prot, residues='all', metadata=True)

def evaluate_steric_clashes(
    ligand_mol: Chem.Mol, receptor_mol: Chem.Mol,
    threshold_overlap: float = 0.4, ifp: Optional[IFP] = None
) -> List[Dict[str, Any]]:
    """Derives clashes from ProLIF's VdWContact interaction (which flags any
    atom pair within van der Waals contact range); overlap_A is computed
    from that contact's exact distance using the same VDW_RADII table as
    before, so downstream severity thresholds (>0.8 Å = severe) are
    unchanged.

    threshold_overlap is a hard cutoff on a receptor treated as rigid: a
    0.41 A overlap is flagged, a 0.39 A overlap is not, with nothing
    graded in between. Docking does not model side-chain relaxation, so a
    slight overlap that a real side chain would rotate away from is
    scored identically to one that could not be resolved at all. This is
    a known limitation of the rigid-receptor approximation generally, not
    specific to this threshold's exact value, and is not addressed here.
    """
    ifp = ifp if ifp is not None else build_ifp(ligand_mol, receptor_mol)
    clashes = []

    for (_, _), interactions in ifp.items():
        for meta in interactions.get('VdWContact', ()):
            l_idx = meta['parent_indices']['ligand'][0]
            p_idx = meta['parent_indices']['protein'][0]
            atom_l = ligand_mol.GetAtomWithIdx(l_idx)
            atom_r = receptor_mol.GetAtomWithIdx(p_idx)
            if atom_l.GetAtomicNum() == 1 or atom_r.GetAtomicNum() == 1:
                continue

            r_vdw_l = get_vdw_radius(atom_l.GetSymbol())
            r_vdw_r = get_vdw_radius(atom_r.GetSymbol())
            dist = meta['distance']
            overlap = (r_vdw_l + r_vdw_r) - dist

            if overlap > threshold_overlap:
                clashes.append({
                    'ligand_atom': f"{atom_l.GetSymbol()}{l_idx}",
                    'receptor_atom': f"{atom_r.GetSymbol()}{p_idx}",
                    'distance_A': round(dist, 2),
                    'overlap_A': round(overlap, 2)
                })

    return sorted(clashes, key=lambda x: x['overlap_A'], reverse=True)

def find_hydrogen_bonds(
    ligand_mol: Chem.Mol, receptor_mol: Chem.Mol, ifp: Optional[IFP] = None
) -> List[Dict[str, Any]]:
    """Ligand<->receptor H-bonds in both directions, from ProLIF's
    HBDonor (ligand donates) / HBAcceptor (ligand accepts) interactions."""
    ifp = ifp if ifp is not None else build_ifp(ligand_mol, receptor_mol)
    hbonds = []

    for (_, pres), interactions in ifp.items():
        for meta in interactions.get('HBDonor', ()):
            hbonds.append({
                'type': 'L_Donor -> R_Acceptor',
                'ligand_atom': f"idx{meta['parent_indices']['ligand'][0]}",
                'receptor_atom': f"{pres}",
                'distance_A': round(meta['distance'], 2),
                'angle_deg': round(meta['DHA_angle'], 1)
            })
        for meta in interactions.get('HBAcceptor', ()):
            hbonds.append({
                'type': 'R_Donor -> L_Acceptor',
                'ligand_atom': f"idx{meta['parent_indices']['ligand'][0]}",
                'receptor_atom': f"{pres}",
                'distance_A': round(meta['distance'], 2),
                'angle_deg': round(meta['DHA_angle'], 1)
            })

    return hbonds

def find_aromatic_contacts(
    ligand_mol: Chem.Mol, receptor_mol: Chem.Mol, ifp: Optional[IFP] = None
) -> List[Dict[str, Any]]:
    """Pi-stacking contacts from ProLIF's FaceToFace/EdgeToFace
    interactions (replaces the old manual ring-centroid/normal geometry)."""
    ifp = ifp if ifp is not None else build_ifp(ligand_mol, receptor_mol)
    aromatic = []

    for (_, pres), interactions in ifp.items():
        for geometry in ('FaceToFace', 'EdgeToFace'):
            for meta in interactions.get(geometry, ()):
                aromatic.append({
                    'distance_A': round(meta['distance'], 2),
                    'geometry': geometry,
                    'receptor_atom': f"{pres}"
                })

    return aromatic


# First-row and related catalytic transition metals posegate.receptor_prep
# retains in a prepared receptor (see receptor_prep.METAL_IONS, duplicated
# here by element symbol rather than imported to keep this a pure geometry
# check with no dependency on how the receptor was prepared).
METAL_ELEMENTS = {'Zn', 'Mg', 'Mn', 'Fe', 'Cu', 'Co', 'Ni', 'Cd', 'Ca'}
# Elements that commonly donate a lone pair to a coordinated metal.
METAL_COORDINATING_ELEMENTS = {'N', 'O', 'S'}
# Typical first-shell metal-donor distance is 2.0-2.3 A; 2.6 A leaves
# margin for docked-pose geometry error without reaching into what would
# clearly be a non-coordinating contact.
METAL_COORDINATION_CUTOFF = 2.6


def find_metal_coordination(
    ligand_mol: Chem.Mol, receptor_mol: Chem.Mol,
    cutoff: float = METAL_COORDINATION_CUTOFF
) -> List[Dict[str, Any]]:
    """Ligand N/O/S atoms within coordination distance of a receptor metal.

    ProLIF has no coordination-bond interaction type, and the metal ions
    posegate.receptor_prep retains are deliberately added without bonds
    (see that module), so coordination is not visible anywhere in the ProLIF
    fingerprint this module otherwise relies on. It is detected here purely
    geometrically instead. For carbonic anhydrase, whose defining
    pharmacophore is a Zn-sulfonamide coordination bond at roughly 2.0 A,
    this is the difference between reporting the actual mechanism and
    reporting nothing more specific than a van der Waals contact.

    This is reporting-only. It is not part of posegate_score: adding a
    sixth feature would require refitting every target's weights, and nothing
    in the five-target validation currently depends on it (the carbonic
    anhydrase result used there rests on the Thr199 hydrogen bond, not on
    zinc coordination).
    """
    metal_atoms = [a for a in receptor_mol.GetAtoms() if a.GetSymbol() in METAL_ELEMENTS]
    if not metal_atoms or ligand_mol.GetNumConformers() == 0 or receptor_mol.GetNumConformers() == 0:
        return []

    rconf = receptor_mol.GetConformer()
    lconf = ligand_mol.GetConformer()
    contacts = []
    for metal in metal_atoms:
        mpos = rconf.GetAtomPosition(metal.GetIdx())
        info = metal.GetPDBResidueInfo()
        metal_label = (f"{info.GetResidueName().strip()}{info.GetResidueNumber()}.{info.GetChainId()}"
                       if info else metal.GetSymbol())
        for atom in ligand_mol.GetAtoms():
            if atom.GetSymbol() not in METAL_COORDINATING_ELEMENTS:
                continue
            dist = lconf.GetAtomPosition(atom.GetIdx()).Distance(mpos)
            if dist <= cutoff:
                contacts.append({
                    'metal': metal_label,
                    'ligand_atom': f"{atom.GetSymbol()}{atom.GetIdx()}",
                    'distance_A': round(dist, 2),
                })
    return sorted(contacts, key=lambda c: c['distance_A'])


def normalize_conserved_residues(
    residues: Union[ConservedResidue, Sequence[ConservedResidue]]
) -> List[ConservedResidue]:
    """Accepts either a single (name, number, chain) tuple or a sequence of
    them, and returns a list. A bare 3-tuple of (str, int, str) is
    ambiguous with a sequence of residues, so it is detected explicitly."""
    if (
        len(residues) == 3
        and isinstance(residues[0], str)
        and isinstance(residues[1], int)
        and isinstance(residues[2], str)
    ):
        return [tuple(residues)]  # type: ignore[list-item]
    return [tuple(r) for r in residues]  # type: ignore[misc]


def find_conserved_hbond(
    ligand_mol: Chem.Mol, receptor_mol: Chem.Mol,
    residues: Union[ConservedResidue, Sequence[ConservedResidue], None] = None,
    mode: str = 'any',
    ifp: Optional[IFP] = None
) -> List[Dict[str, Any]]:
    """Checks for H-bonds to one or more named, conserved receptor residues.

    Args:
        residues: a single (name, number, chain) tuple, or a sequence of
            them. Defaults to BRD4's Asn140. Several targets' defining
            pharmacophores are not a single residue: estrogen receptor
            alpha's charge clamp is Glu353 plus Arg394, and HIV-1
            protease's catalytic dyad is Asp25 on each chain of the
            homodimer.
        mode: 'any' (default) counts the constraint satisfied when at
            least one listed residue is H-bonded; 'all' requires every
            listed residue to be H-bonded, and returns no hbonds unless
            they all are.

    Note that a multi-residue constraint in 'any' mode is easier to
    satisfy than a single-residue one, and in 'all' mode harder. When
    comparing this feature's behaviour across targets whose pharmacophores
    differ in size, that difference in base rate is a property of the
    constraint, not of the poses.
    """
    if mode not in ('any', 'all'):
        raise ValueError(f"mode must be 'any' or 'all', got {mode!r}")

    residue_list = normalize_conserved_residues(
        residues if residues is not None else DEFAULT_CONSERVED_RESIDUES
    )
    ifp = ifp if ifp is not None else build_ifp(ligand_mol, receptor_mol)

    per_residue: Dict[ConservedResidue, List[Dict[str, Any]]] = {r: [] for r in residue_list}

    for (_, pres), interactions in ifp.items():
        for residue in residue_list:
            residue_name, residue_number, chain_id = residue
            if not (pres.name == residue_name and pres.number == residue_number
                    and pres.chain == chain_id):
                continue
            for meta in interactions.get('HBDonor', ()):
                per_residue[residue].append({
                    'type': f'L_Donor -> {residue_name}{residue_number}',
                    'residue': f'{residue_name}{residue_number}.{chain_id}',
                    'distance_A': round(meta['distance'], 2),
                    'angle_deg': round(meta['DHA_angle'], 1)
                })
            for meta in interactions.get('HBAcceptor', ()):
                per_residue[residue].append({
                    'type': f'{residue_name}{residue_number} Donor -> L_Acceptor',
                    'residue': f'{residue_name}{residue_number}.{chain_id}',
                    'distance_A': round(meta['distance'], 2),
                    'angle_deg': round(meta['DHA_angle'], 1)
                })

    if mode == 'all' and not all(per_residue[r] for r in residue_list):
        return []

    return [hb for r in residue_list for hb in per_residue[r]]


# Provisional, BRD4-only weights from scripts/recalibrate_weights.py's fit
# on the 90-compound BRD4 benchmark. See the longer comment inside
# generate_autopsy_report for what these are and are not safe to use for.
POSEGATE_SCORE_WEIGHTS = {
    'hbond_count': 1.047,
    'clash_count': 5.294,
    'aromatic_count': 6.584,
    'conserved_hbond': -2.063,
}


def compute_posegate_score(
    vina_score: float, n_hbonds: int, n_clashes: int, n_aromatic: int, conserved_hit: bool
) -> float:
    """Applies POSEGATE_SCORE_WEIGHTS to raw counts. Pulled out of
    generate_autopsy_report as a pure function so the weights' arithmetic
    (in particular, which terms add vs. subtract) is directly unit-testable
    without needing to build a real docked complex."""
    w = POSEGATE_SCORE_WEIGHTS
    score = vina_score
    score += w['hbond_count'] * n_hbonds
    score += w['clash_count'] * n_clashes
    score += w['aromatic_count'] * n_aromatic
    if conserved_hit:
        score += w['conserved_hbond']
    return score


def generate_autopsy_report(
    ligand_sdf_path: str, receptor_pdb_path: str, vina_score: float,
    conserved_residues: Union[ConservedResidue, Sequence[ConservedResidue], None] = None,
    conserved_mode: str = 'any'
) -> Dict[str, Any]:
    ligand_mol = load_pose_mol(ligand_sdf_path)
    if receptor_pdb_path.endswith('.pkl'):
        # Preferred path for real multi-residue receptors: bonds computed
        # by PDBFixer/OpenMM's Topology and preserved losslessly via
        # pickle, rather than re-guessed from PDB text (see
        # posegate.receptor_prep for why that re-guessing is unreliable).
        receptor_mol = load_receptor_mol(receptor_pdb_path)
    else:
        # proximityBonding=False: RDKit's PDB parser otherwise supplements
        # template-based backbone bonding with a distance-based heuristic
        # that can misfire on tight turns, spuriously bonding two spatially
        # close but unconnected atoms (e.g. residue i's carbonyl C to
        # residue i+1's carbonyl O) and breaking sanitization. Fine for
        # small/synthetic single-residue test receptors; for real
        # multi-residue proteins, prefer the .pkl path above.
        receptor_mol = Chem.MolFromPDBFile(receptor_pdb_path, removeHs=False, proximityBonding=False)

    if ligand_mol is None or receptor_mol is None:
        raise ValueError("Failed to load molecules. Ensure files exist and contain hydrogens.")

    ifp = build_ifp(ligand_mol, receptor_mol)

    report = {
        'vina_score': vina_score,
        'hbonds': find_hydrogen_bonds(ligand_mol, receptor_mol, ifp=ifp),
        'conserved_hbond': find_conserved_hbond(
            ligand_mol, receptor_mol,
            residues=conserved_residues, mode=conserved_mode, ifp=ifp
        ),
        'aromatic': find_aromatic_contacts(ligand_mol, receptor_mol, ifp=ifp),
        # Geometric, not from the ProLIF fingerprint (ProLIF has no
        # coordination-bond type); reporting-only, not part of
        # posegate_score below. See find_metal_coordination.
        'metal_coordination': find_metal_coordination(ligand_mol, receptor_mol),
        'clashes': evaluate_steric_clashes(ligand_mol, receptor_mol, ifp=ifp),
        'decision': 'PENDING',
        'posegate_score': 0.0,
        'explanation': []
    }

    # --- DECISION LOGIC ---
    #
    # PROVISIONAL, BRD4-ONLY DEFAULTS. These per-feature weights come from
    # scripts/recalibrate_weights.py, fit on the 90-compound BRD4 benchmark
    # (30 actives, 60 property-matched decoys) built by the current
    # scripts/fetch_benchmark_dataset.py, predicting active/decoy from
    # (vina_score, hbond_count, conserved_hbond, aromatic_count,
    # clash_count), weights read off relative to the fitted vina_score
    # coefficient so the formula stays in vina_score's own kcal/mol-like
    # units. Cross-validated (out-of-fold) AUC-ROC on that benchmark: 0.619,
    # vs 0.604 for raw vina_score alone.
    #
    # scripts/compare_feature_weights.py, run across BRD4/CDK2/estrogen
    # receptor alpha/HIV-1 protease/carbonic anhydrase, shows that
    # conserved_hbond is the only one of these five features whose sign is
    # bootstrap-stable in the same direction on every target (0.94-1.00
    # sign stability over 200 resamples per target); vina_score, hbond_count
    # and clash_count all reverse sign on at least one target. These BRD4
    # weights are therefore not safe to use unmodified on another target --
    # call scripts/recalibrate_weights.py on that target's own benchmark
    # instead, exactly as the five-target validation did.
    posegate_score = compute_posegate_score(
        vina_score=report['vina_score'],
        n_hbonds=len(report['hbonds']),
        n_clashes=len(report['clashes']),
        n_aromatic=len(report['aromatic']),
        conserved_hit=bool(report['conserved_hbond']),
    )

    # Severe clashes (>0.8 A overlap) are a hard structural-validity gate,
    # not a statistically-fit feature: this is a physically implausible
    # pose regardless of what the calibrated score says about it.
    severe_clashes = [c for c in report['clashes'] if c['overlap_A'] > 0.8]

    report['posegate_score'] = round(posegate_score, 2)

    # Absolute thresholds for standalone single-pose use (no batch to rank
    # against): set to the same BRD4 benchmark's 30th/70th posegate_score
    # percentiles under the weights above, matching rank_batch()'s default
    # 0.3/0.4 fractions so a standalone call and a batch call land on
    # roughly the same operating point. Provisional and BRD4-specific for
    # the same reason the weights are. For screening many candidates at
    # once, prefer rank_batch() (relative ranking) over these fixed cutoffs.
    if severe_clashes:
        report['decision'] = 'REJECT'
        report['explanation'].append(f"Severe steric clash detected (max overlap: {severe_clashes[0]['overlap_A']} Å).")
    elif posegate_score <= -8.66:
        report['decision'] = 'PRIORITIZE'
        report['explanation'].append(f"Strong binding profile (PoseGate Score: {posegate_score}).")
    elif posegate_score <= -7.64:
        report['decision'] = 'REVIEW'
        report['explanation'].append(f"Moderate binding profile (PoseGate Score: {posegate_score}). Requires manual inspection.")
    else:
        report['decision'] = 'REJECT'
        report['explanation'].append(f"Weak binding profile (PoseGate Score: {posegate_score}).")

    return report

def rank_batch(
    reports: List[Dict[str, Any]],
    prioritize_frac: float = 0.3,
    review_frac: float = 0.4
) -> List[Dict[str, Any]]:
    """Re-rank a batch of autopsy reports by percentile of posegate_score,
    for screening many candidates against each other (as in a real virtual
    screen, where you take your top N%) rather than judging a single pose
    against fixed absolute kcal/mol cutoffs, which only make sense relative
    to whatever score distribution a given receptor/scoring setup produces.
    Severe-clash REJECTs from generate_autopsy_report are left untouched;
    everything else is re-decided by rank within this batch. Mutates and
    returns the same report dicts.
    """
    def has_severe_clash(r):
        return any(c['overlap_A'] > 0.8 for c in r['clashes'])

    rankable = [r for r in reports if not has_severe_clash(r)]
    rankable.sort(key=lambda r: r['posegate_score'])

    n = len(rankable)
    n_prioritize = round(n * prioritize_frac)
    n_review = round(n * review_frac)

    for i, r in enumerate(rankable):
        if i < n_prioritize:
            r['decision'] = 'PRIORITIZE'
        elif i < n_prioritize + n_review:
            r['decision'] = 'REVIEW'
        else:
            r['decision'] = 'REJECT'

    return reports
