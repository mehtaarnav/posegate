# posegate/scripts/plip_ensemble_miner.py
"""The same conserved-contact mining posegate.conserved_contacts does, but
sourcing per-structure interactions from PLIP instead of ProLIF, so the two
tools' ensemble-level results can be compared quantitatively rather than by
eye. Mirrors posegate.conserved_contacts.mine_conserved_contacts's output
schema: [{residue, interaction, n_structures, frequency}, ...]."""

import argparse
import json
import os
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict

# PLIP's XML container tag -> a interaction-type label matching, where the
# underlying interaction is comparable, posegate's own ProLIF-derived
# labels (HBDonor/HBAcceptor collapse to one 'HBond' bucket here since
# PLIP's hydrogen_bond entries don't cleanly separate by direction the way
# find_hydrogen_bonds's ligand-role framing does).
PLIP_INTERACTION_TAGS = {
    'hydrogen_bond': 'HBond',
    'hydrophobic_interaction': 'Hydrophobic',
    'pi_stack': 'PiStack',
    'pi_cation_interaction': 'PiCation',
    'salt_bridge': 'SaltBridge',
    'halogen_bond': 'Halogen',
    'water_bridge': 'WaterBridge',
}


def run_plip_full(raw_pdb_path: str, out_dir: str, ligand_resname: str) -> set:
    """Returns {(resnr, chain, interaction_label), ...} for the named
    ligand's binding site only (see compare_plip.py for why: PLIP profiles
    every heteroatom group, including crystallization solvent, as its own
    binding site)."""
    os.makedirs(out_dir, exist_ok=True)
    subprocess.run(
        f"python -m plip.plipcmd -f {raw_pdb_path} -o {out_dir} -x -q",
        shell=True, capture_output=True
    )
    xml_files = [f for f in os.listdir(out_dir) if f.endswith('_report.xml')]
    if not xml_files:
        return set()
    tree = ET.parse(f"{out_dir}/{xml_files[0]}")

    contacts = set()
    for bs in tree.getroot().iter('bindingsite'):
        ident = bs.find('identifiers')
        longname = ident.find('longname') if ident is not None else None
        if longname is None:
            continue
        # PLIP names a binding site by every heterogen group it decided
        # belongs together, e.g. "1SA-ZN" when the ligand sits next to a
        # catalytic metal ion (common for metalloenzymes) -- match if our
        # ligand code is one of the '-'-separated components, not just an
        # exact whole-string match.
        if ligand_resname not in longname.text.split('-'):
            continue
        interactions_el = bs.find('interactions')
        if interactions_el is None:
            continue
        for tag, label in PLIP_INTERACTION_TAGS.items():
            for entry in interactions_el.iter(tag):
                resnr = entry.find('resnr')
                reschain = entry.find('reschain')
                if resnr is None:
                    continue
                chain = reschain.text if reschain is not None else '?'
                contacts.add((int(resnr.text), chain, label))
    return contacts


def mine_plip_ensemble(structures: list, out_dir: str) -> list:
    counts = defaultdict(int)
    n_valid = 0

    for s in structures:
        contacts = run_plip_full(s['pdb_path'], f"{out_dir}/{s['pdb_id']}", s['ligand_resname'])
        if not contacts:
            print(f"  WARNING: no PLIP contacts found for {s['pdb_id']} (ligand {s['ligand_resname']})")
            continue
        n_valid += 1
        for resnr, chain, label in contacts:
            counts[(f"{resnr}.{chain}", label)] += 1

    results = [
        {'residue': residue, 'interaction': interaction, 'n_structures': c, 'frequency': round(c / n_valid, 3)}
        for (residue, interaction), c in counts.items()
    ]
    return sorted(results, key=lambda r: r['frequency'], reverse=True), n_valid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="Raw manifest: [{pdb_id, pdb_path, ligand_resname}, ...]")
    parser.add_argument("--out_dir", default="plip_out")
    parser.add_argument("--out_json", required=True)
    parser.add_argument("--top_n", type=int, default=20)
    args = parser.parse_args()

    with open(args.manifest) as f:
        structures = json.load(f)

    results, n_valid = mine_plip_ensemble(structures, args.out_dir)
    print(f"PLIP-mined {n_valid}/{len(structures)} structures.\n")
    print(f"{'Residue':<14}{'Interaction':<14}{'N':>4}{'Frequency':>12}")
    print('-' * 44)
    for r in results[:args.top_n]:
        print(f"{r['residue']:<14}{r['interaction']:<14}{r['n_structures']:>4}{r['frequency']:>12.2f}")

    with open(args.out_json, 'w') as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
