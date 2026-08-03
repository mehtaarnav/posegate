# posegate/scripts/batch_dock.py
import argparse
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from posegate.docking import prepare_receptor, dock_ligand
from posegate.autopsy import generate_autopsy_report, find_conserved_hbond, rank_batch
from posegate.receptor_prep import (
    prepare_receptor_pickle, load_receptor_mol, write_clean_receptor_pdb
)

def prepare_ligand_sdf(smiles: str, out_path: str):
    """Embeds a 3D conformer, keeping only the largest covalent fragment.

    ChEMBL records many actives as salts, e.g. raloxifene hydrochloride as
    'Cl.O=C(...)'. Embedding those whole produces a multi-fragment molecule
    that OpenBabel writes as PDBQT with a second ROOT, which Vina rejects
    outright ("Unknown or inappropriate tag found in flex residue or
    ligand"), silently dropping the compound from the benchmark. Since the
    counterion is not what binds, the largest fragment is the right thing
    to dock."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    if '.' in smiles:
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        raise ValueError("3D embedding failed")
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


# Per-worker state. The receptor Mol is unpickled once when a worker
# starts rather than once per ligand, since it is the same for every
# ligand and parsing it repeatedly would dominate the per-ligand cost.
_WORKER = {}


def _init_worker(ctx: dict):
    _WORKER.update(ctx)
    _WORKER['receptor_mol'] = load_receptor_mol(ctx['receptor_pkl'])


def _process_ligand(task):
    """Docks and autopsies one ligand. Returns ('ok', report) or
    ('fail', {name, smiles, error}). Runs in a worker process, so it must
    return only picklable data and must not raise."""
    name, smi = task
    try:
        lig_sdf = f"data/{name}.sdf"
        prepare_ligand_sdf(smi, lig_sdf)

        box_size = compute_ligand_box_size(lig_sdf)
        dock_res = dock_ligand(
            _WORKER['receptor_pdbqt'], lig_sdf, _WORKER['center'],
            box_size=box_size, exhaustiveness=_WORKER['exhaustiveness'],
            n_poses=_WORKER['n_poses'], cpu=_WORKER['cpu']
        )

        docked_sdf, selected_score, restraint_met = select_restrained_pose(
            dock_res, _WORKER['receptor_mol'], _WORKER['n_poses'],
            _WORKER['conserved_residues'], _WORKER['conserved_mode']
        )

        report = generate_autopsy_report(
            docked_sdf, _WORKER['receptor_pkl'], selected_score,
            conserved_residues=_WORKER['conserved_residues'],
            conserved_mode=_WORKER['conserved_mode']
        )
        # Stash per-molecule bookkeeping directly on the report dict so
        # it survives the rank_batch() re-ranking pass later.
        report['name'] = name
        report['box_size'] = round(max(box_size), 1)
        report['best_vina_score'] = dock_res['scores'][0]
        report['restraint_met'] = restraint_met
        return ('ok', report)
    except Exception as e:
        return ('fail', {'name': name, 'smiles': smi, 'error': str(e)})


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
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                        help="Ligands to dock concurrently. Vina is pinned to one thread per "
                             "worker, so this trades Vina's internal threading for across-ligand "
                             "parallelism; per-ligand results are unchanged.")
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
    if not os.path.exists(receptor_h_pdb):
        print(f"Preparing heterogen-free hydrogenated receptor: {receptor_h_pdb}")
        write_clean_receptor_pdb(args.receptor_pdb, receptor_h_pdb)
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

    # Ligands are independent, so they are docked concurrently. Vina is
    # pinned to one thread per worker (see dock_ligand's cpu argument) so
    # the workers share the machine instead of each claiming all of it.
    # Every output path is derived from the ligand name, so concurrent
    # workers never write to the same file.
    ctx = {
        'receptor_pdbqt': receptor_pdbqt,
        'receptor_pkl': receptor_pkl,
        'center': args.center,
        'exhaustiveness': args.exhaustiveness,
        'n_poses': args.n_poses,
        'conserved_residues': conserved_residues,
        'conserved_mode': args.conserved_mode,
        'cpu': 1 if args.workers > 1 else 0,
    }
    tasks = [(row['name'], row['smiles']) for _, row in df.iterrows()]
    total = len(tasks)
    print(f"Docking {total} ligands across {args.workers} worker(s)...")

    def checkpoint():
        # Checkpoint raw (pre-ranking) progress so a later crash doesn't
        # lose everything; final decisions still require the full batch.
        pd.DataFrame([
            {'name': r['name'], 'vina_score': r['vina_score'], 'posegate_score': r['posegate_score']}
            for r in reports
        ]).to_csv("results_raw_checkpoint.csv", index=False)

    if args.workers > 1:
        with ProcessPoolExecutor(max_workers=args.workers,
                                 initializer=_init_worker, initargs=(ctx,)) as pool:
            futures = {pool.submit(_process_ligand, t): t[0] for t in tasks}
            for i, fut in enumerate(as_completed(futures), 1):
                status, payload = fut.result()
                if status == 'ok':
                    reports.append(payload)
                else:
                    print(f"  FAILED: {payload['name']} ({payload['error']})")
                    failures.append(payload)
                print(f"[{i}/{total}] {futures[fut]}", flush=True)
                checkpoint()
    else:
        _init_worker(ctx)
        for i, t in enumerate(tasks, 1):
            status, payload = _process_ligand(t)
            if status == 'ok':
                reports.append(payload)
            else:
                print(f"  FAILED: {payload['name']} ({payload['error']})")
                failures.append(payload)
            print(f"[{i}/{total}] {t[0]}", flush=True)
            checkpoint()

    # Restore input order, so a run's output does not depend on the order
    # in which workers happened to finish.
    order = {name: i for i, (name, _) in enumerate(tasks)}
    reports.sort(key=lambda r: order[r['name']])

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
