# posegate/scripts/run_conserved_contact_miner.py
import argparse
import json

from posegate.conserved_contacts import mine_conserved_contacts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="Prepped-structure manifest from prep_ensemble.py")
    parser.add_argument("--top_n", type=int, default=20)
    parser.add_argument("--out_json", default=None)
    args = parser.parse_args()

    with open(args.manifest) as f:
        structures = json.load(f)

    print(f"Mining conserved contacts across {len(structures)} structures: "
          f"{[s['pdb_id'] for s in structures]}")

    results = mine_conserved_contacts(structures)

    print(f"\n{'Residue':<14}{'Interaction':<14}{'N':>4}{'Frequency':>12}")
    print('-' * 44)
    for r in results[:args.top_n]:
        print(f"{r['residue']:<14}{r['interaction']:<14}{r['n_structures']:>4}{r['frequency']:>12.2f}")

    if args.out_json:
        with open(args.out_json, 'w') as f:
            json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
