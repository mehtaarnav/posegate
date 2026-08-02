# posegate/scripts/batch_dock.py
import argparse
import os
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from posegate.docking import prepare_receptor, dock_ligand
from posegate.autopsy import generate_autopsy_report, find_conserved_hbond, rank_batch
from posegate.receptor_prep import prepare_receptor_pickle, load_receptor_mol

def prepare_ligand_sdf(smiles: str, out_path: str):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    Chem.MolToMolFile(mol, out_path)

def compute_ligand_box_size(sdf_path: str, margin: float = 8.0, min_size: float = 16.0):
    """Size the search box to the ligand's own conformer extent plus a
    fixed margin for translational/rotational search room, rather than one
    fixed box for every ligand (a box much smaller than the ligand's own
    extent causes catastrophic clashes; see the JQ1 box-size sweep)."""
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    coords = mol.GetConformer().GetPositions()
    extent = coords.max(axis=0) - coords.min(axis=0)
    return [max(float(e) + margin, min_size) for e in extent]

def select_restrained_pose(dock_res: dict, receptor_mol, n_poses: int,
                            conserved_residues: list, conserved_mode: str = 'any'):
    """Vina's scoring function has no restraint term, so this applies the
    conserved-contact restraint as a pose *selection* filter after the
    fact: among Vina's top n_poses candidates, pick the best-scoring one
    that actually satisfies it, falling back to the single best-scoring
    pose (restraint unmet) if none do. Takes the already-loaded receptor
    Mol (see posegate.receptor_prep) rather than a path, since this is
    called once per ligand and re-parsing from disk each time would be
    wasteful. conserved_residues is a list of (residue_name,
    residue_number, chain_id) -- target-specific (e.g. BRD4's Asn140,
    CDK2's Leu83, or the two residues of estrogen receptor alpha's
    Glu353/Arg394 charge clamp), so it must be supplied by the caller."""
    pose_file = dock_res['pose_file']
    split_prefix = pose_file.replace('.pdbqt', '_p')
    subprocess.run(f"obabel {pose_file} -O {split_prefix}.sdf -m", shell=True, capture_output=True)

    for i, score in enumerate(dock_res['scores']):
        pose_sdf = f"{split_prefix}{i + 1}.sdf"
        if not os.path.exists(pose_sdf):
            continue
        lig = Chem.MolFromMolFile(pose_sdf, removeHs=False)
        if lig is None:
            continue
        if find_conserved_hbond(lig, receptor_mol, residues=conserved_residues,
                                 mode=conserved_mode):
            return pose_sdf, score, True

    fallback_sdf = f"{split_prefix}1.sdf"
    return fallback_sdf, dock_res['scores'][0], False


def parse_conserved_residue(spec: str) -> tuple:
    """Parses a NAME:NUMBER:CHAIN CLI spec, e.g. 'GLU:353:A'."""
    parts = spec.split(':')
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"expected NAME:NUMBER:CHAIN (e.g. GLU:353:A), got {spec!r}")
    name, number, chain = parts
    try:
        return (name.upper(), int(number), chain)
    except ValueError:
        raise argparse.ArgumentTypeError(f"residue number must be an integer in {spec!r}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor_pdb", required=True)
    parser.add_argument("--ligands_csv", required=True, help="CSV with 'name' and 'smiles' columns")
    parser.add_argument("--center", nargs=3, type=float, required=True, help="X Y Z center")
    parser.add_argument("--exhaustiveness", type=int, default=32, help="Vina search exhaustiveness")
    parser.add_argument("--n_poses", type=int, default=9, help="Vina poses per ligand, for restraint-guided selection")
    parser.add_argument("--conserved_residue", nargs='+', type=parse_conserved_residue,
                        default=[('ASN', 140, 'A')], metavar='NAME:NUMBER:CHAIN',
                        help="Conserved-contact residue(s), e.g. ASN:140:A for BRD4, or "
                             "GLU:353:A ARG:394:A for estrogen receptor alpha's charge clamp. "
                             "Supply the residue(s) the miner reported for your target.")
    parser.add_argument("--conserved_mode", choices=['any', 'all'], default='any',
                        help="With several conserved residues, whether a pose must contact any "
                             "of them (default) or all of them")
    args = parser.parse_args()
    conserved_residues = args.conserved_residue

    # Dock against the same heterogen-free, hydrogenated receptor used for
    # autopsy (args.receptor_pdb is the raw crystal file and still contains
    # the original co-crystallized ligand/waters as ATOM records; docking
    # against it silently blocks the real binding site with a "ghost" copy
    # of the native ligand).
    receptor_h_pdb = args.receptor_pdb.replace('.pdb', '_h.pdb')
    receptor_pdbqt = args.receptor_pdb.replace('.pdb', '.pdbqt')
    prepare_receptor(receptor_h_pdb, receptor_pdbqt)

    # Prepared once via posegate.receptor_prep (bonds taken directly from
    # PDBFixer/OpenMM's Topology, not re-guessed from PDB text) and reused
    # for every ligand's autopsy/restraint check below.
    receptor_pkl = args.receptor_pdb.replace('.pdb', '_h.pkl')
    prepare_receptor_pickle(receptor_h_pdb, receptor_pkl)
    receptor_mol = load_receptor_mol(receptor_pkl)

    df = pd.read_csv(args.ligands_csv)
    reports = []
    failures = []

    for _, row in df.iterrows():
        name, smi = row['name'], row['smiles']
        print(f"Processing {name}...")

        try:
            lig_sdf = f"data/{name}.sdf"
            prepare_ligand_sdf(smi, lig_sdf)

            box_size = compute_ligand_box_size(lig_sdf)
            dock_res = dock_ligand(
                receptor_pdbqt, lig_sdf, args.center,
                box_size=box_size, exhaustiveness=args.exhaustiveness, n_poses=args.n_poses
            )

            docked_sdf, selected_score, restraint_met = select_restrained_pose(
                dock_res, receptor_mol, args.n_poses, conserved_residues, args.conserved_mode
            )

            report = generate_autopsy_report(
                docked_sdf, receptor_pkl, selected_score,
                conserved_residues=conserved_residues, conserved_mode=args.conserved_mode
            )
            # Stash per-molecule bookkeeping directly on the report dict so
            # it survives the rank_batch() re-ranking pass below.
            report['name'] = name
            report['box_size'] = round(max(box_size), 1)
            report['best_vina_score'] = dock_res['scores'][0]
            report['restraint_met'] = restraint_met
            reports.append(report)
        except Exception as e:
            print(f"  FAILED: {name} ({e})")
            failures.append({'name': name, 'smiles': smi, 'error': str(e)})

        # Checkpoint raw (pre-ranking) progress so a later crash doesn't
        # lose everything; final decisions still require the full batch.
        pd.DataFrame([
            {'name': r['name'], 'vina_score': r['vina_score'], 'posegate_score': r['posegate_score']}
            for r in reports
        ]).to_csv("results_raw_checkpoint.csv", index=False)

    # Absolute per-molecule thresholds only make sense relative to whatever
    # score distribution this receptor/scoring setup actually produces, so
    # the final PRIORITIZE/REVIEW/REJECT call is made by percentile rank
    # across this batch (severe-clash REJECTs are preserved as-is).
    rank_batch(reports)

    results = [{
        'name': r['name'],
        'box_size': r['box_size'],
        'vina_score': r['vina_score'],
        'best_vina_score': r['best_vina_score'],
        'restraint_met': r['restraint_met'],
        'posegate_score': r['posegate_score'],
        'decision': r['decision'],
        'hbond_count': len(r['hbonds']),
        'conserved_hbond': bool(r['conserved_hbond']),
        'aromatic_count': len(r['aromatic']),
        'clash_count': len(r['clashes'])
    } for r in reports]

    pd.DataFrame(results).to_csv("results.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv("failures.csv", index=False)

    print(f"Docking and Autopsy complete. {len(results)} succeeded, {len(failures)} failed.")
    print("Results saved to results.csv" + (", failures to failures.csv" if failures else ""))

if __name__ == "__main__":
    main()
