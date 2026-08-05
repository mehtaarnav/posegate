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
        lig = Chem.MolFromMolFile(s['ligand_sdf'], removeHs=False)
        receptor_path = s['receptor_pdb']
        if receptor_path.endswith('.pkl'):
            rec = load_receptor_mol(receptor_path)
        else:
            rec = Chem.MolFromPDBFile(receptor_path, removeHs=False, proximityBonding=False)
        if lig is None or rec is None:
            skipped.append(s['pdb_id'])
            continue

        n_valid += 1
        ifp = build_ifp(lig, rec)

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
