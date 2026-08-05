# posegate/scripts/loo_ensemble_size_curve.py
"""Leave-one-out miner accuracy as a function of ensemble size.

Growing the CDK2 mining ensemble from 6 to 19 structures previously made
the miner's output *worse*: ASP145 dropped out of the frequency table
entirely, because 14 of the 19 came from a single fragment-screen series
that only touches the hinge (see scripts/compare_visgremlin.py and the
README's discussion of that comparison). That was diagnosed from one
pair of ensembles, informally. This turns it into a curve: for a range of
ensemble sizes, draw several random subsets of that size from the full
pool, run leave-one-out validation (loo_validate_miner.run_loo) on each,
and report accuracy against size.

If ensemble composition genuinely matters more than count, this curve
should be noisy and non-monotonic rather than a clean increasing curve,
and the variance across random subsets at a fixed size should itself be
informative: high variance at a given size means *which* structures you
have matters more than *how many*.
"""

import argparse
import json
import random
import statistics

# Sibling script, not a package module: importable directly because
# running `python scripts/loo_ensemble_size_curve.py` puts this file's
# own directory (scripts/) on sys.path.
from loo_validate_miner import run_loo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True,
                        help="Prepped-structure manifest covering the full pool of structures "
                             "to subsample from")
    parser.add_argument("--sizes", type=int, nargs='+', default=[6, 10, 14, 18, 22])
    parser.add_argument("--n_subsets", type=int, default=15,
                        help="Random subsets to draw at each size (ignored when size equals "
                             "the full pool, which has only one possible subset)")
    parser.add_argument("--top_k", type=int, default=1,
                        help="Report accuracy at this single top-k value (keeps the sweep to "
                             "one number per subset instead of one per top_k)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out_json", default=None)
    args = parser.parse_args()

    with open(args.manifest) as f:
        pool = json.load(f)
    print(f"Pool: {len(pool)} structures ({args.manifest})")

    rng = random.Random(args.seed)
    results = []
    for size in args.sizes:
        if size > len(pool):
            print(f"size={size}: skipped, exceeds pool of {len(pool)}")
            continue
        # size == len(pool) has exactly one possible subset (the whole
        # pool, order doesn't matter to run_loo), so repeating the draw
        # produces n_subsets copies of the same result, not independent
        # trials. Run it once and report a single point with no spread.
        n_trials = 1 if size == len(pool) else args.n_subsets
        accs = []
        seen_subsets = set()
        for trial in range(n_trials):
            subset = rng.sample(pool, size)
            key = frozenset(s['pdb_id'] for s in subset)
            if key in seen_subsets and size < len(pool):
                # Small pool, large size: distinct subsets are scarce
                # enough that resampling can collide. Redraw once rather
                # than silently double-counting the same subset.
                subset = rng.sample(pool, size)
                key = frozenset(s['pdb_id'] for s in subset)
            seen_subsets.add(key)
            fold_results = run_loo(subset, [args.top_k])
            if not fold_results:
                continue
            hits = sum(r[f'top{args.top_k}_hit'] for r in fold_results)
            acc = hits / len(fold_results)
            accs.append(acc)
        if not accs:
            print(f"size={size}: no valid trials")
            continue
        mean = statistics.mean(accs)
        spread = statistics.pstdev(accs) if len(accs) > 1 else 0.0
        print(f"size={size:2d}  mean_top{args.top_k}_accuracy={mean:.2f}  "
              f"stdev={spread:.2f}  trials={[round(a, 2) for a in accs]}")
        results.append({'size': size, 'mean_accuracy': mean, 'stdev': spread, 'trials': accs})

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()
