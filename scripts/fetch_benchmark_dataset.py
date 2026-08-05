# posegate/scripts/fetch_benchmark_dataset.py
"""Fetches real actives from ChEMBL for any target, plus DUD-E-style
property-matched decoys, parameterized by target_chembl_id. Used to build
the benchmark for each of the five validated targets, so
scripts/recalibrate_weights.py's fitted weights can be compared across
them (see the cross-target feature-weight comparison in
scripts/compare_feature_weights.py)."""

import argparse
import os

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

def fetch_candidate_pool(pool_size=20000, page_size=1000):
    """Pull a large, target-agnostic pool of ChEMBL molecules as decoy candidates.

    The ChEMBL API caps a single response at 1000 records, so this pages
    through until pool_size is reached. The pool has to be large: decoys
    are accepted only if they match an active on five properties at once,
    and a target whose actives sit at the edge of chemical space will
    exhaust a small pool and fall back on whatever poorly-matched
    candidates remain. Estrogen receptor alpha is the motivating case,
    with actives averaging logP 5.76, which a 1000-molecule pool could not
    match without skewing the benchmark.
    """
    candidates = []
    for offset in range(0, pool_size, page_size):
        url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json"
               f"?limit={min(page_size, pool_size - offset)}&offset={offset}")
        try:
            molecules = requests.get(url, timeout=60).json().get('molecules', [])
        except Exception as e:
            print(f"  pool fetch stopped at offset {offset}: {e}")
            break
        if not molecules:
            break
        for mol in molecules:
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
        'hbd': Lipinski.NumHDonors(mol),
        'hba': Lipinski.NumHAcceptors(mol),
        'fp': AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
    }

def property_matched_decoys(
    actives, pool_smiles,
    mw_tol=0.25, logp_tol=1.5, rotb_tol=2, hbd_tol=1, hba_tol=2,
    sim_cutoff=0.35, per_active=2
):
    """DUD-E-style decoy selection: property-matched to an active, but
    structurally dissimilar (low Tanimoto) so it isn't just a near-analog.

    Hydrogen-bond donor and acceptor counts are matched alongside
    molecular weight, logP and rotatable bonds. Without them, actives and
    decoys can differ systematically in how many hydrogen bonds they are
    capable of forming at all, which would make posegate's hbond_count
    feature reflect the benchmark's construction rather than anything
    about binding. Since the transferability analysis turns on whether
    that feature changes sign between targets, leaving donors and
    acceptors unmatched would confound exactly the comparison being made.

    Dissimilarity is enforced against *every* active, not just the active
    currently being matched. Checking only the current active lets a
    molecule that happens to be a near-analog of some other active enter
    the decoy set, where it may well be a genuine binder.
    """
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

    all_active_fps = [d['fp'] for d in active_descs]

    # Precompute each candidate's similarity to its most similar active, so
    # the dissimilarity filter costs one pass rather than one per active.
    for c in pool_descs:
        c['max_sim'] = max(DataStructs.BulkTanimotoSimilarity(c['fp'], all_active_fps))

    def property_distance(candidate, active):
        """Normalized distance across the matched properties, so no single
        property dominates through its units. Selecting the nearest
        candidates rather than the first ones inside the tolerance window
        matters when a target's actives sit at the edge of the pool's
        distribution: every acceptable candidate then lies on the same side
        of the active, and taking them in pool order biases the whole decoy
        set that way. Estrogen receptor alpha showed this on logP, with
        decoys landing 0.8 log units below the actives while still inside
        the window."""
        return (
            abs(candidate['mw'] - active['mw']) / max(active['mw'] * mw_tol, 1e-9)
            + abs(candidate['logp'] - active['logp']) / logp_tol
            + abs(candidate['rotb'] - active['rotb']) / max(rotb_tol, 1e-9)
            + abs(candidate['hbd'] - active['hbd']) / max(hbd_tol, 1e-9)
            + abs(candidate['hba'] - active['hba']) / max(hba_tol, 1e-9)
        )

    used_smiles = set()
    decoys = []
    for a_desc in active_descs:
        eligible = [
            c for c in pool_descs
            if c['smiles'] not in used_smiles
            and abs(c['mw'] - a_desc['mw']) <= a_desc['mw'] * mw_tol
            and abs(c['logp'] - a_desc['logp']) <= logp_tol
            and abs(c['rotb'] - a_desc['rotb']) <= rotb_tol
            and abs(c['hbd'] - a_desc['hbd']) <= hbd_tol
            and abs(c['hba'] - a_desc['hba']) <= hba_tol
            and c['max_sim'] <= sim_cutoff
        ]
        eligible.sort(key=lambda c: property_distance(c, a_desc))

        for c in eligible[:per_active]:
            decoys.append({'name': f"Decoy_{len(decoys)}", 'smiles': c['smiles'], 'ic50_nM': None})
            used_smiles.add(c['smiles'])
    return decoys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target_chembl_id", required=True)
    parser.add_argument("--limit", type=int, default=50, help="Raw activities to fetch before filtering")
    parser.add_argument("--per_active", type=int, default=2, help="Property-matched decoys per active")
    parser.add_argument("--pool_size", type=int, default=20000,
                        help="Decoy candidate pool size. Larger pools give better property "
                             "matching for targets whose actives are unusual; check the result "
                             "with scripts/check_decoy_quality.py")
    parser.add_argument("--pool_cache", default=None,
                        help="File to cache the decoy candidate pool in (one SMILES per line). "
                             "Reused across targets so every benchmark draws decoys from the "
                             "same pool.")
    parser.add_argument("--out_csv", required=True)
    args = parser.parse_args()

    print(f"Fetching actives for {args.target_chembl_id} from ChEMBL...")
    actives = fetch_chembl_actives(args.target_chembl_id, limit=args.limit)
    print(f"Found {len(actives)} actives.")

    # The decoy pool is target-agnostic, so cache it: building benchmarks
    # for several targets should not re-download the same molecules each
    # time, and reusing one pool also keeps the decoy sets comparable
    # across targets, which matters when the fitted weights are then
    # compared between them.
    pool = []
    if args.pool_cache and os.path.exists(args.pool_cache):
        with open(args.pool_cache) as f:
            pool = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(pool)} pool candidates from {args.pool_cache}")

    # Only refetch on a cache miss. Comparing against pool_size would
    # refetch every time, since the API returns fewer usable molecules
    # than requested (records without structures are dropped), so a cache
    # of 19878 would never satisfy a request for 20000.
    if not pool:
        print("Fetching decoy candidate pool from ChEMBL...")
        pool = fetch_candidate_pool(pool_size=args.pool_size)
        if args.pool_cache:
            with open(args.pool_cache, 'w') as f:
                f.write('\n'.join(pool))
            print(f"Cached pool to {args.pool_cache}")
    print(f"Pool size: {len(pool)}")

    print("Property-matching decoys (MW/logP/rotatable bonds/donors/acceptors matched, "
          "Tanimoto-dissimilar to every active)...")
    decoys = property_matched_decoys(actives, pool, per_active=args.per_active)
    print(f"Selected {len(decoys)} property-matched decoys.")

    df = pd.DataFrame(actives + decoys)
    df.to_csv(args.out_csv, index=False)
    print(f"Saved {args.out_csv}")

if __name__ == "__main__":
    main()
