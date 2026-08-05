# posegate/posegate/conserved_contacts.py
"""PDB-ensemble conserved-contact mining.

Given a set of PDB structures for the same target, each co-crystallized
with a different ligand, aggregates ProLIF interaction fingerprints across
the ensemble to surface which receptor contacts are conserved (present
across many structures/ligands) versus incidental to one specific ligand.

This generalizes the previously hand-picked, hardcoded constraint in
posegate.autopsy.find_conserved_hbond (BRD4's Asn140, chosen from reading
the literature) into a data-driven, target-agnostic pipeline: point it at
any target's set of ligand-bound PDB structures and it surfaces the
equivalent conserved-contact residues automatically, from the structures
themselves rather than a hardcoded residue name/number.

Output rows are per (residue, interaction type), and VdWContact is a
superset of the specific types: any hydrogen-bonded, hydrophobic or
aromatic contact is necessarily also within van der Waals range, so a
residue can appear as both e.g. 'VdWContact 1.00' and 'HBDonor 0.80' with
nothing in either row indicating the second is a more specific
description of (part of) the first. Callers that want only the specific
interaction types should filter 'VdWContact' out themselves (see
scripts/run_conserved_contact_miner.py's --exclude_vdw, and
scripts/compare_visgremlin.py's POSEGATE_SPECIFIC, which already does
this for its own comparison).
"""

import math
from collections import defaultdict
from typing import Any, Dict, List, Tuple

from rdkit import Chem

from posegate.autopsy import build_ifp
from posegate.receptor_prep import load_receptor_mol


def _load_structure(structure: Dict[str, str]):
    """Loads a structure's ligand and receptor, returning (None, None) on
    any failure rather than letting a missing file's OSError propagate
    differently from RDKit's own None-on-malformed-input convention.
    Chem.MolFromMolFile raises OSError for a file that does not exist but
    returns None for one that exists and fails to parse; every caller
    here already treats a None ligand/receptor as 'skip this structure',
    so an uncaught OSError previously crashed the whole ensemble instead
    of skipping the one bad structure -- including every leave-one-out
    fold whose *training* set happened to include it, not just a fold
    that held it out directly."""
    try:
        lig = Chem.MolFromMolFile(structure['ligand_sdf'], removeHs=False)
    except OSError:
        lig = None
    try:
        receptor_path = structure['receptor_pdb']
        if receptor_path.endswith('.pkl'):
            rec = load_receptor_mol(receptor_path)
        else:
            rec = Chem.MolFromPDBFile(receptor_path, removeHs=False, proximityBonding=False)
    except OSError:
        rec = None
    return lig, rec


def wilson_interval(count: int, n: int, z: float = 1.96) -> Tuple[float, float]:
    """95%-default Wilson score interval for a binomial proportion.

    A raw frequency treats 3/6 and 11/22 as the same number, 50%, but
    they carry very different statistical weight: the first is consistent
    with anywhere from about 19% to 81% under a 95% interval, the second
    with about 31% to 69%. mine_conserved_contacts' frequency threshold
    convention (commonly read at 0.5) has no such distinction built in, so
    this interval is reported alongside every frequency rather than
    silently treating a 6-structure and a 22-structure ensemble as
    equally conclusive at the same cutoff. Wilson's interval is used
    rather than the normal approximation because it stays inside [0, 1]
    and remains reasonable at small n and at frequencies near 0 or 1,
    both of which are common here (a residue seen in 1 of 6 structures).
    """
    if n == 0:
        return (0.0, 0.0)
    p = count / n
    denom = 1 + z * z / n
    center = p + z * z / (2 * n)
    adj = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (center - adj) / denom), min(1.0, (center + adj) / denom))


def mine_conserved_contacts(structures: List[Dict[str, str]]) -> List[Dict[str, Any]]:
    """Mines conserved receptor contacts across a PDB ensemble.

    Args:
        structures: a list of dicts, one per ensemble member, each with:
            - 'pdb_id': identifier, for reporting only
            - 'ligand_sdf': path to that structure's bound ligand (with Hs)
            - 'receptor_pdb': path to that structure's receptor, prepared
              via posegate.receptor_prep.prepare_receptor_pickle (.pkl,
              preferred) or a plain heterogen-free hydrogenated PDB (.pdb,
              fallback; see posegate.autopsy.generate_autopsy_report for
              why .pkl is preferred for real multi-residue receptors)

    Returns:
        A list of {'residue', 'interaction', 'n_structures', 'n_ensemble',
        'frequency', 'ci95'} dicts, sorted by descending frequency (fraction
        of the ensemble in which that (residue, interaction type) pair
        occurs at least once). A residue interacting with every ligand's
        own chemically distinct scaffold, via the same interaction type,
        is a conserved contact; one hit in a single structure is very
        likely specific to that one ligand rather than a general
        pharmacophore feature of the pocket.

        'frequency' is a point estimate, not a statistically established
        value: a threshold like 0.5 means something different for a
        6-structure ensemble than for a 22-structure one, and no
        correction is applied for scoring many residues from the same
        ensemble at once. 'ci95' (a Wilson score interval; see
        wilson_interval) makes that uncertainty explicit rather than
        implicit in the ensemble size alone.
    """
    counts: Dict[tuple, int] = defaultdict(int)
    n_valid = 0
    skipped = []

    for s in structures:
        lig, rec = _load_structure(s)
        if lig is None or rec is None:
            skipped.append(s['pdb_id'])
            continue

        # ProLIF's interaction detection can fail on a structure that
        # loaded successfully -- e.g. VdWContact raises ValueError for an
        # element its chosen radii table has no entry for, which real
        # structures do contain (old-style heavy-atom phasing derivatives
        # such as mercury or platinum, retained in the deposited
        # coordinates even though they play no role in binding). Uncaught,
        # this crashes mining for the whole ensemble over one structure;
        # skipped and reported like a load failure instead.
        try:
            ifp = build_ifp(lig, rec)
        except Exception as e:
            skipped.append(f"{s['pdb_id']} (interaction detection failed: {e})")
            continue

        n_valid += 1

        # Count each (residue, interaction type) at most once per
        # structure, so a residue with many contacting atoms in one
        # structure doesn't outweigh a residue seen across many structures.
        seen_this_structure = set()
        for (_, pres), interactions in ifp.items():
            for iname in interactions:
                key = (str(pres), iname)
                seen_this_structure.add(key)
        for key in seen_this_structure:
            counts[key] += 1

    if n_valid == 0:
        raise ValueError("No structures could be loaded; check ligand_sdf/receptor_pdb paths.")
    if skipped:
        print(f"Skipped {len(skipped)}/{len(structures)} structure(s) that failed to parse "
              f"(RDKit rejected the ligand or receptor): {skipped}")

    results = [
        {
            'residue': residue,
            'interaction': interaction,
            'n_structures': count,
            'n_ensemble': n_valid,
            'frequency': round(count / n_valid, 3),
            'ci95': tuple(round(x, 3) for x in wilson_interval(count, n_valid)),
        }
        for (residue, interaction), count in counts.items()
    ]
    return sorted(results, key=lambda r: r['frequency'], reverse=True)


# Interaction types specific enough to count as a real contact when
# checking a held-out structure below, as opposed to VdWContact (see this
# module's docstring on why that one is excluded by default elsewhere).
SPECIFIC_INTERACTIONS = {
    'HBDonor', 'HBAcceptor', 'Hydrophobic', 'FaceToFace', 'EdgeToFace', 'PiStacking'
}


def _top_k_predicted_residues(mined_rows: List[Dict[str, Any]], k: int,
                               exclude_vdw: bool = True) -> List[str]:
    """The miner's top-k distinct residues by descending frequency.
    mined_rows is one (residue, interaction) row per line, already sorted;
    this collapses to distinct residues, since a held-out check is
    per-residue, not per-interaction-type."""
    seen: List[str] = []
    for row in mined_rows:
        if exclude_vdw and row['interaction'] == 'VdWContact':
            continue
        residue = row['residue'].upper()
        if residue not in seen:
            seen.append(residue)
        if len(seen) >= k:
            break
    return seen


def _structure_contact_residues(structure: Dict[str, str]):
    """Residues a single structure's own ligand specifically contacts,
    computed directly rather than via mining. Returns None if the
    structure fails to load, or if it loads but ProLIF's own interaction
    detection fails on it (see mine_conserved_contacts for why that is a
    real, separate failure mode from a load failure -- a bound heavy atom
    ProLIF's radii table has no entry for is the case actually seen)."""
    lig, rec = _load_structure(structure)
    if lig is None or rec is None:
        return None

    try:
        ifp = build_ifp(lig, rec)
    except Exception:
        return None
    residues = set()
    for (_, pres), interactions in ifp.items():
        if SPECIFIC_INTERACTIONS.intersection(interactions):
            residues.add(str(pres).upper())
    return residues


# Ensemble-size reliability thresholds, derived from an actual measured
# curve, not a guess. A 22-structure CDK2 ensemble was leave-one-out
# validated at sizes 6/10/14/18/22, 15 random subsets per size (except
# 22, the full pool, which has only one possible subset):
#
#   size   mean accuracy   stdev
#   6      0.59            0.23   (individual trials ranged 0.0-1.0)
#   10     0.71            0.15
#   14     0.67            0.05
#   18     0.69            0.06
#   22     0.68            --
#
# Accuracy stops improving past ~10 structures, but the variance is what
# actually matters for trusting a single run: at 6 structures, two draws
# of the same size landed at 0% and 100% accuracy. By 14+ the spread
# tightens to roughly 0.6-0.8. This threshold is derived from one target
# (CDK2); it is the best evidence available, not a guarantee it transfers
# exactly to every target's geometry.
RELIABILITY_THRESHOLDS = (
    (10, 'low', "fewer than 10 structures: in the measured curve this size regime had a "
                "3x higher spread than larger ensembles (stdev 0.15-0.23 vs 0.05-0.06), "
                "including a same-size draw that scored 0% and another that scored 100%. "
                "Treat this accuracy as a rough signal, not a precise estimate."),
    (14, 'moderate', "10-13 structures: variance was still meaningfully elevated in the "
                     "measured curve (stdev ~0.15) compared to 14+ (~0.05). Usable, but "
                     "expect this number to move if you add or swap a few structures."),
)


def ensemble_reliability(n_usable: int) -> Dict[str, Any]:
    """Classifies an ensemble size against RELIABILITY_THRESHOLDS.
    Returns {'tier', 'note'}; tier is 'low', 'moderate', or 'high'."""
    for threshold, tier, note in RELIABILITY_THRESHOLDS:
        if n_usable < threshold:
            return {'tier': tier, 'note': note}
    return {'tier': 'high', 'note': "14+ structures: in the measured curve this is the "
                                     "range where accuracy stabilizes (stdev 0.05-0.06)."}


def leave_one_out_validate(
    structures: List[Dict[str, str]], top_k: Tuple[int, ...] = (1, 3, 5)
) -> Dict[str, Any]:
    """Self-validates the miner on its own input ensemble. Needs no ground
    truth beyond the ensemble the caller already supplied.

    Every validation of the miner up to this point compared its output,
    mined from an ensemble, against a literature pharmacophore known
    before that ensemble was built -- a claim of agreement with prior
    knowledge, and a circular one, since the answer was known in advance.
    This asks a harder question that does not require knowing the answer
    in advance: for each structure in the ensemble in turn, mine the
    remaining N-1, take the miner's top-k predicted residues, and check
    whether the held-out structure's own ligand actually contacts them.
    The held-out structure was never seen by the fold that predicted it.

    Because this only consumes the ensemble the caller already has, it
    runs unchanged on a target with no known literature pharmacophore --
    the exact situation the miner exists for -- and produces a confidence
    number for *this* ensemble specifically, rather than an appeal to
    validation performed on some other, previously studied target.
    run_conserved_contact_miner.py runs this automatically on every
    invocation for that reason (see --skip_self_validation there to opt
    out on a very large ensemble, where it costs one extra full mining
    pass per structure).

    A prediction counts as a hit if the predicted residue shows any
    specific (non-van-der-Waals) interaction with the held-out ligand --
    the same residue, not necessarily the same interaction type, since
    the interaction type recorded on the training folds is what informed
    the prediction and need not match exactly what the held-out ligand's
    own chemistry produces.

    Returns a dict with 'folds' (one entry per structure, with its
    held-out contacts, the miner's top-k predictions, and a hit/miss flag
    per k) and 'accuracy' (per k: hit count, usable fold count, and
    accuracy), plus 'n_ensemble' and 'n_usable' for structures that failed
    to load and were skipped.
    """
    folds = []
    for i, held in enumerate(structures):
        train = structures[:i] + structures[i + 1:]
        actual = _structure_contact_residues(held)
        if actual is None:
            folds.append({'pdb_id': held['pdb_id'], 'skipped': True})
            continue

        mined = mine_conserved_contacts(train)
        fold: Dict[str, Any] = {
            'pdb_id': held['pdb_id'],
            'skipped': False,
            'held_out_residues': sorted(actual),
        }
        for k in top_k:
            predicted = _top_k_predicted_residues(mined, k)
            fold[f'top{k}_predicted'] = predicted
            fold[f'top{k}_hit'] = any(p in actual for p in predicted)
        folds.append(fold)

    usable = [f for f in folds if not f['skipped']]
    accuracy = {}
    for k in top_k:
        hits = sum(f[f'top{k}_hit'] for f in usable)
        n = len(usable)
        accuracy[k] = {
            'hits': hits,
            'n_folds': n,
            'accuracy': round(hits / n, 3) if n else None,
        }

    return {
        'folds': folds,
        'accuracy': accuracy,
        'n_ensemble': len(structures),
        'n_usable': len(usable),
        'reliability': ensemble_reliability(len(usable)),
    }
