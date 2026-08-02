# posegate/scripts/compare_plip.py
"""Cross-validates posegate's ProLIF-based H-bond detection against PLIP
(an independent, widely-used interaction-detection tool) on the same real
structures, across 5 diverse protein families. Agreement is compared at
the residue level (which receptor residues are flagged as H-bonding to
the ligand), since exact atom-level bookkeeping differs between tools."""

import argparse
import json
import re
import subprocess
import xml.etree.ElementTree as ET

from rdkit import Chem

from posegate.autopsy import find_hydrogen_bonds
from posegate.receptor_prep import load_receptor_mol


def run_plip(raw_pdb_path: str, out_dir: str, ligand_resname: str) -> set:
    """Residues H-bonding to the *named* ligand only.

    PLIP profiles every heteroatom group in the file as its own binding
    site, including crystallization additives (DMSO, ethylene glycol,
    glycerol, ions). Collecting hydrogen_bond elements from the whole
    report therefore mixes solvent contacts in with the real ligand's --
    e.g. 3MXF yields 5 sites (JQ1 + DMS + 3x EDO), and only one of them
    is the ligand of interest. Filter to the matching binding site.
    """
    subprocess.run(
        f"python -m plip.plipcmd -f {raw_pdb_path} -o {out_dir} -x -q",
        shell=True, capture_output=True
    )
    xml_files = [f for f in __import__('os').listdir(out_dir) if f.endswith('_report.xml')]
    if not xml_files:
        return set()
    tree = ET.parse(f"{out_dir}/{xml_files[0]}")

    residues = set()
    matched = False
    for bs in tree.getroot().iter('bindingsite'):
        ident = bs.find('identifiers')
        longname = ident.find('longname') if ident is not None else None
        # PLIP names a binding site by every heterogen group it decided
        # belongs together (e.g. "1SA-ZN" next to a catalytic metal ion),
        # so match on component membership, not exact equality.
        if longname is None or ligand_resname not in longname.text.split('-'):
            continue
        matched = True
        for hb in bs.iter('hydrogen_bond'):
            residues.add(int(hb.find('resnr').text))
    if not matched:
        print(f"  WARNING: no PLIP binding site found for ligand {ligand_resname}")
    return residues


def run_posegate(ligand_sdf: str, receptor_pkl: str) -> set:
    lig = Chem.MolFromMolFile(ligand_sdf, removeHs=False)
    rec = load_receptor_mol(receptor_pkl)
    hbonds = find_hydrogen_bonds(lig, rec)
    residues = set()
    for hb in hbonds:
        m = re.search(r'(\d+)', hb['receptor_atom'])
        if m:
            residues.add(int(m.group(1)))
    return residues


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="JSON: [{pdb_id, raw_pdb, ligand_sdf, receptor_pkl}, ...]")
    parser.add_argument("--out_dir", default="plip_out")
    args = parser.parse_args()

    with open(args.manifest) as f:
        targets = json.load(f)

    rows = []
    for t in targets:
        pdb_id = t['pdb_id']
        print(f"Comparing {pdb_id}...")
        plip_out = f"{args.out_dir}/{pdb_id}"
        __import__('os').makedirs(plip_out, exist_ok=True)

        plip_residues = run_plip(t['raw_pdb'], plip_out, t['ligand_resname'])
        posegate_residues = run_posegate(t['ligand_sdf'], t['receptor_pkl'])

        both = plip_residues & posegate_residues
        either = plip_residues | posegate_residues
        jaccard = len(both) / len(either) if either else None

        rows.append({
            'pdb_id': pdb_id,
            'plip_residues': sorted(plip_residues),
            'posegate_residues': sorted(posegate_residues),
            'agree': sorted(both),
            'plip_only': sorted(plip_residues - posegate_residues),
            'posegate_only': sorted(posegate_residues - plip_residues),
            'jaccard': round(jaccard, 3) if jaccard is not None else None,
        })

    print()
    for r in rows:
        print(f"{r['pdb_id']}:")
        print(f"  PLIP residues:     {r['plip_residues']}")
        print(f"  posegate residues: {r['posegate_residues']}")
        print(f"  agree: {r['agree']}  |  PLIP-only: {r['plip_only']}  |  posegate-only: {r['posegate_only']}")
        print(f"  Jaccard overlap: {r['jaccard']}")
        print()

    with open(f"{args.out_dir}/comparison.json", 'w') as f:
        json.dump(rows, f, indent=2)


if __name__ == "__main__":
    main()
