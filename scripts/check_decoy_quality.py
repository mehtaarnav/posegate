# posegate/scripts/check_decoy_quality.py
"""Quality control for a benchmark built by fetch_benchmark_dataset.py.

Decoy quality is a confound for any cross-target comparison of fitted
feature weights. If decoys are property-matched to actives on one target
but not on another, a feature can appear to change its relationship to
activity when what actually changed is the decoy set. This reports, per
benchmark, whether the DUD-E-style matching actually held.

Two failure modes matter:

  1. Property skew. Decoys are supposed to match actives on molecular
     weight, logP and rotatable-bond count, so that those bulk properties
     cannot themselves separate the classes. Standardized mean difference
     (Cohen's d) is reported per property; |d| >= 0.5 is flagged, since a
     medium effect in a bulk property means the benchmark is partly
     separable without any structural reasoning at all.
  2. Insufficient dissimilarity. Decoys are supposed to be structurally
     dissimilar to actives (Tanimoto below the selection cutoff), so they
     are not near-analogs that might genuinely bind. The maximum
     similarity of any decoy to any active is reported.

Hydrogen-bond donor and acceptor counts are reported for a further
reason: posegate's transferability analysis turns on whether the
hbond_count feature reverses sign between targets, and that feature's
behaviour is directly downstream of how well matched donors and
acceptors are between the classes.
"""

import argparse

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski

RDLogger.DisableLog('rdApp.*')

MATCHED_PROPERTIES = ['mw', 'logp', 'rotb']
REPORTED_PROPERTIES = MATCHED_PROPERTIES + ['hbd', 'hba']


def describe(smiles: str):
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


def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    pooled = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
    if pooled == 0:
        return 0.0
    return float((a.mean() - b.mean()) / pooled)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark_csv", required=True)
    parser.add_argument("--sim_cutoff", type=float, default=0.35,
                        help="Tanimoto cutoff used when the decoys were selected")
    parser.add_argument("--d_threshold", type=float, default=0.5,
                        help="Flag a matched property whose |Cohen's d| reaches this")
    args = parser.parse_args()

    df = pd.read_csv(args.benchmark_csv)
    df['is_active'] = df['name'].astype(str).str.startswith('Active')

    described, labels = [], []
    for _, row in df.iterrows():
        d = describe(row['smiles'])
        if d is not None:
            described.append(d)
            labels.append(bool(row['is_active']))
    labels = np.array(labels)

    n_act, n_dec = int(labels.sum()), int((~labels).sum())
    print(f"{args.benchmark_csv}: {len(described)} parsed, {n_act} actives, {n_dec} decoys")
    if n_act == 0 or n_dec == 0:
        raise SystemExit("Benchmark has no actives or no decoys.")
    print(f"active:decoy ratio 1:{n_dec / n_act:.1f}\n")

    print(f"{'property':<8}{'actives':>10}{'decoys':>10}{'Cohen d':>10}   status")
    flagged = []
    for prop in REPORTED_PROPERTIES:
        vals = np.array([d[prop] for d in described], dtype=float)
        a, b = vals[labels], vals[~labels]
        d = cohens_d(a, b)
        if prop in MATCHED_PROPERTIES:
            ok = abs(d) < args.d_threshold
            status = 'matched' if ok else 'SKEWED'
            if not ok:
                flagged.append(prop)
        else:
            status = 'reported only'
        print(f"{prop:<8}{a.mean():>10.2f}{b.mean():>10.2f}{d:>10.2f}   {status}")

    active_fps = [d['fp'] for d, is_a in zip(described, labels) if is_a]
    max_sims = [
        max(DataStructs.BulkTanimotoSimilarity(d['fp'], active_fps))
        for d, is_a in zip(described, labels) if not is_a
    ]
    over = sum(s > args.sim_cutoff for s in max_sims)
    print(f"\ndecoy max-Tanimoto to any active: mean {np.mean(max_sims):.3f}, "
          f"max {np.max(max_sims):.3f}")
    print(f"decoys above the {args.sim_cutoff} selection cutoff: {over}")

    duplicates = len(df[~df['is_active']]) - df[~df['is_active']]['smiles'].nunique()
    print(f"duplicate decoy structures: {duplicates}")

    print()
    if flagged:
        print(f"FAIL: property-matching did not hold for {', '.join(flagged)}. These bulk "
              f"properties partly separate actives from decoys on their own, so fitted "
              f"feature weights from this benchmark are not comparable against a target "
              f"whose matching did hold.")
    elif over:
        print(f"WARN: {over} decoy(s) exceed the similarity cutoff and may be genuine binders.")
    else:
        print("PASS: matched properties are balanced and decoys are structurally dissimilar.")


if __name__ == "__main__":
    main()
