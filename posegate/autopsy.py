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
    unchanged."""
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

def generate_autopsy_report(
    ligand_sdf_path: str, receptor_pdb_path: str, vina_score: float,
    conserved_residues: Union[ConservedResidue, Sequence[ConservedResidue], None] = None,
    conserved_mode: str = 'any'
) -> Dict[str, Any]:
    ligand_mol = Chem.MolFromMolFile(ligand_sdf_path, removeHs=False)
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
        'clashes': evaluate_steric_clashes(ligand_mol, receptor_mol, ifp=ifp),
        'decision': 'PENDING',
        'posegate_score': 0.0,
        'explanation': []
    }

    # --- DECISION LOGIC ---
    #
    # These per-feature weights come from scripts/recalibrate_weights.py:
    # an L2-regularized logistic regression fit on the 65-compound BRD4
    # benchmark (22 actives, 43 property-matched decoys), predicting
    # active/decoy from (vina_score, hbond_count, conserved_hbond,
    # aromatic_count, clash_count), with weights read off relative to the
    # fitted vina_score coefficient so the whole formula stays in
    # vina_score's own kcal/mol-like units. Cross-validated (out-of-fold)
    # AUC-ROC on that benchmark: 0.618, vs 0.542 for raw vina_score alone
    # — a real, non-circular improvement, but calibrated on one target and
    # 65 compounds; treat these specific numbers as provisional until
    # validated on more targets, not as universal constants.
    #
    # Generic hbond_count got a positive (penalizing) weight: in this
    # benchmark, decoys tended to form more generic/incidental H-bonds
    # than actives, while the specific conserved_hbond contact remained a
    # strong reward. An earlier version of this comment attributed that to
    # the decoys being property-matched on donor/acceptor counts, which
    # they were not: the decoy selection in scripts/fetch_brd4_dataset.py
    # matched only molecular weight, logP and rotatable bonds, so actives
    # and decoys were free to differ in how many H-bonds they could form
    # at all. scripts/fetch_benchmark_dataset.py now also matches donor
    # and acceptor counts, so these weights should be refitted from a
    # benchmark rebuilt with that script before they are relied on.
    #
    # aromatic_count was regularized to exactly zero in every fit, which
    # was previously read as the feature not being useful. It was in fact
    # never measured: prepared receptors had no aromatic rings at all
    # until posegate.receptor_prep began assigning side-chain aromaticity,
    # so ProLIF could not report pi-stacking and the count was constant at
    # zero. Nothing is known about whether this feature helps; it needs a
    # refit against receptors prepared by the current code.
    posegate_score = report['vina_score']
    posegate_score += 2.706 * len(report['hbonds'])
    posegate_score += 0.813 * len(report['clashes'])
    if report['conserved_hbond']:
        posegate_score -= 3.104

    # Severe clashes (>0.8 A overlap) are a hard structural-validity gate,
    # not a statistically-fit feature: this is a physically implausible
    # pose regardless of what the calibrated score says about it.
    severe_clashes = [c for c in report['clashes'] if c['overlap_A'] > 0.8]

    report['posegate_score'] = round(posegate_score, 2)

    # Absolute thresholds for standalone single-pose use (no batch to rank
    # against): set to the same benchmark's 30th/70th posegate_score
    # percentiles, matching rank_batch()'s default 0.3/0.4 fractions so a
    # standalone call and a batch call land on roughly the same operating
    # point. For screening many candidates at once, prefer rank_batch()
    # (relative ranking) over these fixed cutoffs.
    if severe_clashes:
        report['decision'] = 'REJECT'
        report['explanation'].append(f"Severe steric clash detected (max overlap: {severe_clashes[0]['overlap_A']} Å).")
    elif posegate_score <= -8.5:
        report['decision'] = 'PRIORITIZE'
        report['explanation'].append(f"Strong binding profile (PoseGate Score: {posegate_score}).")
    elif posegate_score <= -6.7:
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
