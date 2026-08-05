# posegate/scripts/run_conserved_contact_miner.py
import argparse
import json

from posegate.conserved_contacts import mine_conserved_contacts


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
    parser.add_argument("--out_json", default=None)
    args = parser.parse_args()

    with open(args.manifest) as f:
        structures = json.load(f)

    print(f"Mining conserved contacts across {len(structures)} structures: "
          f"{[s['pdb_id'] for s in structures]}")

    results = mine_conserved_contacts(structures)
    if args.exclude_vdw:
        results = [r for r in results if r['interaction'] != 'VdWContact']

    print(f"\n{'Residue':<14}{'Interaction':<14}{'N':>4}{'Frequency':>12}  95% CI")
    print('-' * 60)
    for r in results[:args.top_n]:
        lo, hi = r['ci95']
        print(f"{r['residue']:<14}{r['interaction']:<14}{r['n_structures']:>4}{r['frequency']:>12.2f}"
              f"  [{lo:.2f}, {hi:.2f}]")

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
