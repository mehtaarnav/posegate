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
"""

from collections import defaultdict
from typing import Any, Dict, List

from rdkit import Chem

from posegate.autopsy import build_ifp
from posegate.receptor_prep import load_receptor_mol


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
        A list of {'residue', 'interaction', 'n_structures', 'frequency'}
        dicts, sorted by descending frequency (fraction of the ensemble in
        which that (residue, interaction type) pair occurs at least once).
        A residue interacting with every ligand's own chemically distinct
        scaffold, via the same interaction type, is a conserved contact;
        one hit in a single structure is very likely specific to that one
        ligand rather than a general pharmacophore feature of the pocket.
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
            'frequency': round(count / n_valid, 3)
        }
        for (residue, interaction), count in counts.items()
    ]
    return sorted(results, key=lambda r: r['frequency'], reverse=True)
