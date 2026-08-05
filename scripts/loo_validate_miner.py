# posegate/scripts/loo_validate_miner.py
"""CLI for posegate.conserved_contacts.leave_one_out_validate.

The leave-one-out logic itself now lives in the library
(posegate.conserved_contacts.leave_one_out_validate), not here: it is
also run automatically by run_conserved_contact_miner.py on every mining
call, and having two copies of the same validation logic -- one used in
production, one only in this standalone script -- is exactly the kind of
drift that produces a silently stale research script. This file is now
only the command-line entry point and the per-fold progress printing.
"""

import argparse
import json

from posegate.conserved_contacts import leave_one_out_validate


def run_loo(structures, top_k_values):
    """Runs leave-one-out validation, printing per-fold progress, and
    returns the flat list of fold dicts (skipped folds excluded) that
    this script and scripts/loo_ensemble_size_curve.py both consume."""
    result = leave_one_out_validate(structures, top_k=tuple(top_k_values))
    out = []
    for fold in result['folds']:
        if fold.get('skipped'):
            print(f"  {fold['pdb_id']}: SKIPPED (structure failed to load)")
            continue
        summary = "  ".join(
            f"top{k}={'HIT' if fold[f'top{k}_hit'] else 'miss'}" for k in top_k_values
        )
        print(f"  {fold['pdb_id']}: {summary}")
        out.append(fold)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True,
                        help="Prepped-structure manifest from prep_ensemble.py")
    parser.add_argument("--top_k", type=int, nargs='+', default=[1, 3, 5],
                        help="Report hit rate for the miner's top-k predicted residues, "
                             "at each of these k values")
    parser.add_argument("--out_json", default=None)
    args = parser.parse_args()

    with open(args.manifest) as f:
        structures = json.load(f)

    print(f"Leave-one-out validation over {len(structures)} structures "
          f"({args.manifest})")
    results = run_loo(structures, args.top_k)

    n = len(results)
    print(f"\n{n}/{len(structures)} folds completed")
    for k in args.top_k:
        hits = sum(r[f'top{k}_hit'] for r in results)
        print(f"  top-{k} accuracy: {hits / n:.2f}  ({hits}/{n})")

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nWrote {args.out_json}")


if __name__ == "__main__":
    main()
