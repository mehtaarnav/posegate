# posegate/scripts/fetch_benchmark_dataset.py
"""Generalized version of fetch_brd4_dataset.py: fetches real actives from
ChEMBL for any target + DUD-E-style property-matched decoys, parameterized
by target_chembl_id rather than hardcoded to BRD4. Used for the
transferability study (scripts/recalibrate_weights.py's fitted weights
were calibrated only on BRD4; this builds equivalent benchmarks for other
targets to test whether those weights generalize)."""

import argparse
import requests
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski

def fetch_chembl_actives(target_chembl_id, limit=100, ic50_cutoff_nm=1000):
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={target_chembl_id}&limit={limit}&standard_type=IC50"
    response = requests.get(url).json()

    actives = []
    for activity in response['activities']:
        # Only accept exact, nM-normalized values; ChEMBL mixes in
        # inequality relations ("<"/">") and non-nM units that would
        # otherwise be misread as literal potency.
        if (
            activity.get('standard_relation') == '='
            and activity.get('standard_units') == 'nM'
            and activity.get('standard_value')
            and float(activity['standard_value']) < ic50_cutoff_nm
        ):
            smiles = activity.get('canonical_smiles')
            if smiles:
                actives.append({'name': f"Active_{len(actives)}", 'smiles': smiles, 'ic50_nM': activity['standard_value']})
    return actives

def fetch_candidate_pool(pool_size=1000):
    """Pull a large, target-agnostic pool of ChEMBL molecules as decoy candidates."""
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit={pool_size}"
    response = requests.get(url).json()

    candidates = []
    for mol in response['molecules']:
        smiles = (mol.get('molecule_structures') or {}).get('canonical_smiles')
        if smiles:
            candidates.append(smiles)
    return candidates

def compute_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        'mw': Descriptors.MolWt(mol),
        'logp': Crippen.MolLogP(mol),
        'rotb': Lipinski.NumRotatableBonds(mol),
        'fp': AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
    }

def property_matched_decoys(
    actives, pool_smiles,
    mw_tol=0.25, logp_tol=1.5, rotb_tol=2,
    sim_cutoff=0.35, per_active=2
):
    """DUD-E-style decoy selection: property-matched to an active, but
    structurally dissimilar (low Tanimoto) so it isn't just a near-analog."""
    active_smiles = {a['smiles'] for a in actives}
    active_descs = []
    for a in actives:
        d = compute_descriptors(a['smiles'])
        if d:
            active_descs.append(d)

    pool_descs = []
    for smi in pool_smiles:
        if smi in active_smiles:
            continue
        d = compute_descriptors(smi)
        if d:
            d['smiles'] = smi
            pool_descs.append(d)

    used_smiles = set()
    decoys = []
    for a_desc in active_descs:
        matches = 0
        for c in pool_descs:
            if c['smiles'] in used_smiles:
                continue
            if abs(c['mw'] - a_desc['mw']) > a_desc['mw'] * mw_tol:
                continue
            if abs(c['logp'] - a_desc['logp']) > logp_tol:
                continue
            if abs(c['rotb'] - a_desc['rotb']) > rotb_tol:
                continue
            sim = DataStructs.TanimotoSimilarity(a_desc['fp'], c['fp'])
            if sim > sim_cutoff:
                continue

            decoys.append({'name': f"Decoy_{len(decoys)}", 'smiles': c['smiles'], 'ic50_nM': None})
            used_smiles.add(c['smiles'])
            matches += 1
            if matches >= per_active:
                break
    return decoys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target_chembl_id", required=True)
    parser.add_argument("--limit", type=int, default=50, help="Raw activities to fetch before filtering")
    parser.add_argument("--per_active", type=int, default=2, help="Property-matched decoys per active")
    parser.add_argument("--out_csv", required=True)
    args = parser.parse_args()

    print(f"Fetching actives for {args.target_chembl_id} from ChEMBL...")
    actives = fetch_chembl_actives(args.target_chembl_id, limit=args.limit)
    print(f"Found {len(actives)} actives.")

    print("Fetching decoy candidate pool from ChEMBL...")
    pool = fetch_candidate_pool(pool_size=1000)
    print(f"Pool size: {len(pool)}")

    print("Property-matching decoys (MW/logP/rotatable-bonds matched, Tanimoto-dissimilar)...")
    decoys = property_matched_decoys(actives, pool, per_active=args.per_active)
    print(f"Selected {len(decoys)} property-matched decoys.")

    df = pd.DataFrame(actives + decoys)
    df.to_csv(args.out_csv, index=False)
    print(f"Saved {args.out_csv}")

if __name__ == "__main__":
    main()
