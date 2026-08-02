# posegate/scripts/run_autopsy.py
import argparse
from posegate.autopsy import generate_autopsy_report

def print_report(report, name):
    print(f"\n{'='*50}")
    print(f"POSEGATE AUTOPSY: {name}")
    print(f"{'='*50}")
    print(f"Vina Score:  {report['vina_score']} kcal/mol")
    print(f"Decision:    {report['decision']}")
    print(f"{'-'*50}")
    
    if report['hbonds']:
        print(f"H-BONDS ({len(report['hbonds'])}):")
        for h in report['hbonds'][:3]:
            print(f"  {h['type']} | Dist: {h['distance_A']}Å | Ang: {h['angle_deg']}°")
            
    if report['aromatic']:
        print(f"AROMATIC ({len(report['aromatic'])}):")
        for a in report['aromatic'][:3]:
            print(f"  {a['geometry']} | Dist: {a['distance_A']}Å")
            
    if report['clashes']:
        print(f"CLASHES ({len(report['clashes'])}):")
        for c in report['clashes'][:3]:
            print(f"  {c['ligand_atom']} <-> {c['receptor_atom']} | Overlap: {c['overlap_A']}Å")
            
    print(f"{'-'*50}")
    print("EXPLANATION:")
    for e in report['explanation']:
        print(f"  - {e}")
    print(f"{'='*50}\n")

def main():
    parser = argparse.ArgumentParser(description="PoseGate Autopsy CLI")
    parser.add_argument("--receptor", required=True, help="Path to receptor PDB (with H)")
    parser.add_argument("--ligand", required=True, help="Path to ligand SDF (with H)")
    parser.add_argument("--score", type=float, default=0.0, help="Vina score")
    parser.add_argument("--name", default="Unknown", help="Ligand name")
    args = parser.parse_args()

    report = generate_autopsy_report(args.ligand, args.receptor, args.score)
    print_report(report, args.name)

if __name__ == "__main__":
    main()