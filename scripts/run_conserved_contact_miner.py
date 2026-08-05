# posegate/scripts/run_conserved_contact_miner.py
import argparse
import json

from posegate.conserved_contacts import mine_conserved_contacts, leave_one_out_validate


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="Prepped-structure manifest from prep_ensemble.py")
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--exclude_vdw", action="store_true",
                        help="Drop VdWContact rows from the printed/saved output. VdWContact "
                             "fires on any atom pair within van der Waals radii and is a "
                             "superset of the specific interaction types (HBDonor, HBAcceptor, "
                             "Hydrophobic, FaceToFace, EdgeToFace): a hydrogen-bonded residue "
                             "will also show up as a VdWContact, at its own separate frequency, "
                             "and nothing in the row itself says the two are the same contact. "
                             "Use this to see only the specific interaction types.")
    parser.add_argument("--skip_self_validation", action="store_true",
                        help="Skip the automatic leave-one-out self-validation pass. It costs "
                             "one extra full mining run per structure (N total, versus 1), which "
                             "matters on a large ensemble; skip it if you already trust this "
                             "target's result or are iterating quickly.")
    parser.add_argument("--out_json", default=None,
                        help="Writes the mined (residue, interaction) table, unchanged in "
                             "shape from before self-validation existed: a bare list, so "
                             "existing consumers (compare_visgremlin.py, compare_miners.py) "
                             "keep working. Self-validation goes to --self_validation_json "
                             "instead, not folded in here.")
    parser.add_argument("--self_validation_json", default=None,
                        help="Writes the leave-one-out self-validation result (folds + "
                             "accuracy) to this separate file. No effect with "
                             "--skip_self_validation.")
    args = parser.parse_args()

    with open(args.manifest) as f:
        structures = json.load(f)

    print(f"Mining conserved contacts across {len(structures)} structures: "
          f"{[s['pdb_id'] for s in structures]}")

    results = mine_conserved_contacts(structures)
    reported = [r for r in results if not (args.exclude_vdw and r['interaction'] == 'VdWContact')]

    print(f"\n{'Residue':<14}{'Interaction':<14}{'N':>4}{'Frequency':>12}  95% CI")
    print('-' * 60)
    for r in reported[:args.top_n]:
        lo, hi = r['ci95']
        print(f"{r['residue']:<14}{r['interaction']:<14}{r['n_structures']:>4}{r['frequency']:>12.2f}"
              f"  [{lo:.2f}, {hi:.2f}]")

    self_validation = None
    if not args.skip_self_validation:
        # Runs on this same ensemble, with no ground truth beyond it, so
        # it works identically on a target with no known literature
        # pharmacophore -- the answer to "should I trust this?" that a
        # frequency table alone cannot give.
        print(f"\nSelf-validating: leave-one-out over these same {len(structures)} structures...")
        self_validation = leave_one_out_validate(structures)
        print(f"\n{'top-k':<8}{'hits/folds':<14}accuracy")
        print('-' * 32)
        for k, acc in self_validation['accuracy'].items():
            n = acc['n_folds']
            hits = acc['hits']
            val = f"{acc['accuracy']:.2f}" if acc['accuracy'] is not None else "n/a"
            print(f"top-{k:<4}{f'{hits}/{n}':<14}{val}")
        if self_validation['n_usable'] < self_validation['n_ensemble']:
            skipped = self_validation['n_ensemble'] - self_validation['n_usable']
            print(f"({skipped} structure(s) skipped: failed to load)")
        reliability = self_validation['reliability']
        print(f"\nEnsemble size reliability: {reliability['tier'].upper()}")
        print(f"  {reliability['note']}")

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)
    if args.self_validation_json and self_validation is not None:
        with open(args.self_validation_json, 'w') as f:
            json.dump(self_validation, f, indent=2)


if __name__ == "__main__":
    main()
