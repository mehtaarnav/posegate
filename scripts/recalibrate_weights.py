# posegate/scripts/recalibrate_weights.py
"""Recalibrates posegate_score's feature weights against ProLIF's actual
output distribution, instead of the hand-picked constants in
generate_autopsy_report (tuned against an earlier, less complete
hand-rolled interaction detector and never revisited after the ProLIF
migration).

Fits an L2-regularized logistic regression predicting active/decoy from
(vina_score, hbond_count, conserved_hbond, aromatic_count, clash_count),
reports honest out-of-fold (cross-validated) AUC so the result isn't
circular/overfit to this same 65-compound set, and converts the fitted
coefficients into an additive posegate_score formula in the same units
(kcal/mol-like) as vina_score, for a direct comparison against the old
hand-picked weights and against raw vina_score alone.
"""

import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict

FEATURES = ['vina_score', 'hbond_count', 'conserved_hbond', 'aromatic_count', 'clash_count']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_csv", default="results.csv")
    parser.add_argument("--n_folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    df = pd.read_csv(args.results_csv)
    df['label'] = df['name'].apply(lambda n: 1 if str(n).startswith('Active') else 0)
    df['conserved_hbond'] = df['conserved_hbond'].astype(int)

    X = df[FEATURES].to_numpy(dtype=float)
    y = df['label'].to_numpy()

    print(f"N={len(df)}, actives={y.sum()}, decoys={len(y) - y.sum()}, features={FEATURES}")

    # Baseline: raw vina_score alone.
    auc_vina = roc_auc_score(y, -df['vina_score'])
    print(f"\nBaseline AUC (raw vina_score alone):        {auc_vina:.3f}")

    # Fit L2-regularized logistic regression, choosing regularization
    # strength by nested cross-validation (LogisticRegressionCV), then get
    # an honest out-of-fold AUC via cross_val_predict rather than scoring
    # the same data the model was fit on.
    cv = StratifiedKFold(n_splits=args.n_folds, shuffle=True, random_state=args.seed)
    model = LogisticRegressionCV(cv=cv, penalty='l2', scoring='roc_auc', max_iter=5000)

    oof_scores = cross_val_predict(model, X, y, cv=cv, method='decision_function')
    auc_oof = roc_auc_score(y, oof_scores)
    print(f"Cross-validated AUC (recalibrated, out-of-fold): {auc_oof:.3f}")

    # Fit on the full data to report the final weights (used for
    # deployment; the number above, not this fit's in-sample AUC, is the
    # honest estimate of how well it generalizes).
    model.fit(X, y)
    coefs = dict(zip(FEATURES, model.coef_[0]))
    print(f"\nFitted coefficients (higher magnitude = more discriminating): {coefs}")

    # Convert to an additive posegate_score formula in vina_score's own
    # units: divide every other coefficient by the vina_score coefficient
    # so the formula reads as "vina_score plus/minus so many kcal/mol
    # worth of credit per feature", matching the existing architecture.
    c_vina = coefs['vina_score']
    if c_vina == 0:
        print("WARNING: vina_score coefficient is ~0; cannot express as an additive kcal/mol-scale formula.")
        return

    # decision_function z = sum(c_i * x_i) has higher z = more active-like.
    # posegate_score should have *lower* = more active-like (matches
    # vina_score's own convention), so posegate_score ~ -z. Dividing every
    # coefficient by c_vina (not -c_vina) gives exactly that: since
    # vina_score's own coefficient c_vina is expected negative (more
    # negative vina_score -> more active-like -> positive contribution to
    # z), dividing by c_vina rescales the vina term to coefficient +1 and
    # simultaneously flips every other term's sign correctly, so each
    # w below is already the correct signed additive weight.
    print("\nRecalibrated posegate_score formula: posegate_score = vina_score + sum(w_feat * feature)")
    for feat in FEATURES[1:]:
        w = coefs[feat] / c_vina
        sign_note = " (reward, lowers score)" if w < 0 else " (penalty, raises score)"
        print(f"  {feat}: {w:+.3f} per unit{sign_note}")

if __name__ == "__main__":
    main()
