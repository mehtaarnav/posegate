# posegate/scripts/compare_feature_weights.py
"""Cross-target comparison of fitted pose-feature weights.

The transferability question is not whether the fitted score reaches a
good AUC on any one target -- with ~50-220 compounds per target the AUC
confidence intervals are wide enough that no pairwise difference between
targets is significant, so a table of per-target AUCs would not support a
conclusion. The answerable question is categorical: does a given feature
keep the same *direction* of association with activity across targets, or
does it reverse?

Two failure modes make a naive sign comparison misleading, and both are
handled here:

  1. Near-zero coefficients. A feature regularized to near zero has an
     essentially arbitrary sign, so reporting that it "flipped" between
     two targets manufactures a finding out of noise. Sign stability is
     therefore estimated by bootstrap: the fraction of resamples in which
     the coefficient keeps its sign. A feature whose sign is not stable
     within a target cannot be said to agree or disagree across targets,
     and is reported as indeterminate rather than counted either way.
  2. Scale differences between features. Coefficients are fitted on
     standardized features so that magnitudes are comparable across
     features measured in different units (kcal/mol vs counts).

Feeds on the per-target results CSVs written by batch_dock.py.
"""

import argparse
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

FEATURES = ['vina_score', 'hbond_count', 'conserved_hbond', 'aromatic_count', 'clash_count']


def load(results_csv: str):
    df = pd.read_csv(results_csv)
    df['label'] = df['name'].astype(str).str.startswith('Active').astype(int)
    df['conserved_hbond'] = df['conserved_hbond'].astype(int)
    X = df[FEATURES].to_numpy(dtype=float)
    y = df['label'].to_numpy()
    return X, y


def fit_once(X, y, seed):
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    model = LogisticRegressionCV(cv=cv, penalty='l2', scoring='roc_auc', max_iter=5000)
    model.fit(X, y)
    return model.coef_[0]


def analyse(label, results_csv, n_boot, seed):
    X_raw, y = load(results_csv)
    X = StandardScaler().fit_transform(X_raw)

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    model = LogisticRegressionCV(cv=cv, penalty='l2', scoring='roc_auc', max_iter=5000)
    oof = cross_val_predict(model, X, y, cv=cv, method='decision_function')
    auc_fit = roc_auc_score(y, oof)
    # Raw vina_score baseline: lower (more negative) is more active-like.
    auc_vina = roc_auc_score(y, -X_raw[:, 0])

    coefs = fit_once(X, y, seed)

    # Bootstrap for sign stability, stratified so class balance is kept.
    rng = np.random.default_rng(seed)
    idx_pos = np.flatnonzero(y == 1)
    idx_neg = np.flatnonzero(y == 0)
    boot = np.zeros((n_boot, len(FEATURES)))
    for b in range(n_boot):
        take = np.concatenate([
            rng.choice(idx_pos, size=len(idx_pos), replace=True),
            rng.choice(idx_neg, size=len(idx_neg), replace=True),
        ])
        try:
            boot[b] = fit_once(X[take], y[take], seed)
        except Exception:
            boot[b] = np.nan

    stability = {}
    for i, f in enumerate(FEATURES):
        col = boot[:, i]
        col = col[~np.isnan(col)]
        if len(col) == 0:
            stability[f] = 0.5
            continue
        same = np.mean(np.sign(col) == np.sign(coefs[i])) if coefs[i] != 0 else 0.5
        stability[f] = float(same)

    return {
        'label': label, 'n': len(y), 'actives': int(y.sum()),
        'auc_vina': auc_vina, 'auc_fit': auc_fit,
        'coefs': dict(zip(FEATURES, coefs)), 'stability': stability,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", nargs='+', required=True,
                        metavar='LABEL=PATH',
                        help="Per-target results CSVs, e.g. brd4=data/brd4_results.csv")
    parser.add_argument("--n_boot", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--stability_threshold", type=float, default=0.90,
                        help="Bootstrap sign agreement below which a coefficient's direction "
                             "is treated as indeterminate rather than positive or negative")
    args = parser.parse_args()

    results = []
    for spec in args.results:
        label, path = spec.split('=', 1)
        if not os.path.exists(path):
            print(f"skipping {label}: {path} not found")
            continue
        results.append(analyse(label, path, args.n_boot, args.seed))

    if not results:
        raise SystemExit("No results to compare.")

    print(f"\n{'target':<10}{'N':>6}{'actives':>9}{'raw vina':>10}{'fitted':>9}")
    print('-' * 44)
    for r in results:
        print(f"{r['label']:<10}{r['n']:>6}{r['actives']:>9}{r['auc_vina']:>10.3f}{r['auc_fit']:>9.3f}")

    print("\nStandardized coefficients (sign = direction of association with activity)")
    print("An asterisk marks a coefficient whose sign is not stable under bootstrap;")
    print(f"its direction is indeterminate at this sample size (threshold {args.stability_threshold:.2f}).\n")
    header = f"{'feature':<18}" + ''.join(f"{r['label']:>16}" for r in results)
    print(header)
    print('-' * len(header))

    verdicts = {}
    for f in FEATURES:
        row = f"{f:<18}"
        signs = []
        for r in results:
            c = r['coefs'][f]
            stable = r['stability'][f] >= args.stability_threshold
            row += f"{c:>+14.3f}{'' if stable else '*':<2}"
            signs.append(np.sign(c) if stable else 0)
        print(row)

        known = [s for s in signs if s != 0]
        n_stable, n_total = len(known), len(results)
        support = f"{n_stable}/{n_total} targets stable"
        if n_stable < 2:
            verdicts[f] = f'indeterminate ({support}, too few to compare)'
        elif len(set(known)) == 1:
            verdicts[f] = f'consistent ({support})'
        else:
            verdicts[f] = f'REVERSES between targets ({support})'

    print("\nVerdict per feature:")
    for f in FEATURES:
        print(f"  {f:<18} {verdicts[f]}")
    print("\nA verdict rests only on the targets whose sign was stable. 'consistent' on 2 of 5 "
          "\nis far weaker evidence than 'consistent' on 5 of 5, so the support count is given.")

    print("\nBootstrap sign stability:")
    hdr = f"{'feature':<18}" + ''.join(f"{r['label']:>16}" for r in results)
    print(hdr)
    for f in FEATURES:
        print(f"{f:<18}" + ''.join(f"{r['stability'][f]:>16.2f}" for r in results))


if __name__ == "__main__":
    main()
