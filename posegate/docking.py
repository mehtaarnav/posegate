# posegate/posegate/docking.py
import subprocess
from pathlib import Path
from vina import Vina

def prepare_receptor(pdb_path: str, pdbqt_path: str):
    """Uses OpenBabel to convert PDB to PDBQT and add hydrogens."""
    cmd = f"obabel {pdb_path} -O {pdbqt_path} -xr -h"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

def prepare_ligand(ligand_sdf: str, ligand_pdbqt: str):
    """Uses OpenBabel to convert an SDF ligand to PDBQT (Vina requires PDBQT)."""
    cmd = f"obabel {ligand_sdf} -O {ligand_pdbqt}"
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

def dock_ligand(
    receptor_pdbqt: str,
    ligand_sdf: str,
    center: list,
    box_size: list = [20, 20, 20],
    exhaustiveness: int = 8,
    n_poses: int = 1
) -> dict:
    """Runs AutoDock Vina and returns pose score(s).

    Vina's scoring function has no restraint term, so this returns the
    top n_poses candidates (score + multi-model PDBQT) rather than a
    single answer, so a caller can apply a restraint (e.g. a required
    H-bond) as a pose *selection* filter after the fact.
    """
    ligand_pdbqt = ligand_sdf.replace('.sdf', '.pdbqt')
    prepare_ligand(ligand_sdf, ligand_pdbqt)

    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

    scores = [float(e[0]) for e in v.energies(n_poses=n_poses)]

    out_pdbqt = ligand_sdf.replace('.sdf', '_docked.pdbqt')
    v.write_poses(out_pdbqt, n_poses=n_poses, overwrite=True)

    return {'score': scores[0], 'scores': scores, 'pose_file': out_pdbqt}
