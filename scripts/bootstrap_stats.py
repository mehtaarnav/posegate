# posegate/scripts/bootstrap_stats.py
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from rdkit.ML.Scoring.Scoring import CalcBEDROC, CalcEnrichment

def compute_metrics(df: pd.DataFrame, alpha: float = 20.0, ef_fraction: float = 0.01):
    """Point-estimate AUC-ROC, EF at ef_fraction, and BEDROC(alpha) for one
    (possibly resampled) active/decoy set. Lower posegate_score = more
    active-like, so it's ranked ascending / negated for AUC."""
    ranked = df.sort_values('posegate_score', ascending=True)
    scores_col = [[1.0, int(label)] for label in ranked['label']]
    auc = roc_auc_score(ranked['label'], -ranked['posegate_score'])
    ef = CalcEnrichment(scores_col, 1, [ef_fraction])[0]
    bedroc = CalcBEDROC(scores_col, 1, alpha=alpha)
    return auc, ef, bedroc

def stratified_bootstrap(df: pd.DataFrame, n_boot: int, alpha: float, ef_fraction: float, seed: int):
    """Resamples actives and decoys separately (each with replacement, at
    their own original counts) rather than resampling the pooled set, so
    every resample has the same class balance as the real data and the
    resulting CI reflects genuine metric uncertainty rather than variance
    from occasionally drawing a resample with a different active:decoy
    ratio than the actual dataset."""
    rng = np.random.default_rng(seed)
    actives = df[df['label'] == 1]
    decoys = df[df['label'] == 0]

    aucs, efs, bedrocs = [], [], []
    for _ in range(n_boot):
        a_sample = actives.iloc[rng.integers(0, len(actives), size=len(actives))]
        d_sample = decoys.iloc[rng.integers(0, len(decoys), size=len(decoys))]
        sample = pd.concat([a_sample, d_sample], ignore_index=True)
        try:
            auc, ef, bedroc = compute_metrics(sample, alpha=alpha, ef_fraction=ef_fraction)
        except Exception:
            continue
        aucs.append(auc)
        efs.append(ef)
        bedrocs.append(bedroc)

    return np.array(aucs), np.array(efs), np.array(bedrocs)

def ci(vals: np.ndarray):
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_csv", default="results.csv")
    parser.add_argument("--n_boot", type=int, default=5000)
    parser.add_argument("--alpha", type=float, default=20.0, help="BEDROC alpha (early recognition weighting)")
    parser.add_argument("--ef_fraction", type=float, default=0.01, help="EF cutoff fraction, e.g. 0.01 for EF1%%")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = pd.read_csv(args.results_csv)
    df['label'] = df['name'].apply(lambda n: 1 if str(n).startswith('Active') else 0)
    n_actives = int(df['label'].sum())
    n_decoys = len(df) - n_actives
    print(f"N={len(df)}, actives={n_actives}, decoys={n_decoys}")

    top_n = max(1, round(len(df) * args.ef_fraction))
    print(f"Note: top {args.ef_fraction * 100:.0f}% of {len(df)} compounds = {top_n} compound(s); "
          f"EF/BEDROC are low-precision at this N.")

    auc_point, ef_point, bedroc_point = compute_metrics(df, alpha=args.alpha, ef_fraction=args.ef_fraction)
    print(f"Point estimates: AUC={auc_point:.3f}  EF{args.ef_fraction*100:.0f}%={ef_point:.2f}  BEDROC(alpha={args.alpha:.0f})={bedroc_point:.3f}")

    aucs, efs, bedrocs = stratified_bootstrap(df, args.n_boot, args.alpha, args.ef_fraction, args.seed)
    print(f"Valid stratified bootstrap resamples: {len(aucs)}/{args.n_boot}")

    auc_ci, ef_ci, bedroc_ci = ci(aucs), ci(efs), ci(bedrocs)
    print(f"AUC-ROC:              {auc_point:.3f}  95% CI [{auc_ci[0]:.3f}, {auc_ci[1]:.3f}]")
    print(f"EF{args.ef_fraction*100:.0f}%:                 {ef_point:.2f}  95% CI [{ef_ci[0]:.2f}, {ef_ci[1]:.2f}]")
    print(f"BEDROC(alpha={args.alpha:.0f}):        {bedroc_point:.3f}  95% CI [{bedroc_ci[0]:.3f}, {bedroc_ci[1]:.3f}]")

if __name__ == "__main__":
    main()
