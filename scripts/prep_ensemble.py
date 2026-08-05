# posegate/scripts/prep_ensemble.py
"""Prepares a PDB-ensemble for conserved-contact mining: for each raw PDB
file + its bound ligand's residue name, extracts a hydrogenated ligand
SDF and a receptor pickle (see posegate.receptor_prep for why a pickled
RDKit Mol, built from PDBFixer/OpenMM's own Topology bonds, is used
instead of a plain PDB file for the receptor)."""

import argparse
import json
import subprocess

from posegate.receptor_prep import prepare_receptor_pickle
from posegate.residue_mapping import build_residue_map, remap_residue_numbers


def extract_single_ligand_instance(pdb_path: str, resname: str, out_path: str) -> tuple:
    """Extracts exactly one instance of the ligand (the asymmetric unit
    often contains multiple copies of the same complex; grabbing every
    HETATM matching the resname merges them into one nonsensical molecule).
    Returns (n_atoms, chain_id) for the instance that was kept, using the
    first (chain, residue number) encountered as that single instance."""
    key = None
    n = 0
    chain_id = None
    with open(pdb_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('HETATM') and line[17:20].strip() == resname:
                this_key = (line[21], line[22:26])
                if key is None:
                    key = this_key
                    chain_id = line[21]
                if this_key == key:
                    fout.write(line)
                    n += 1
    return n, chain_id


def chains_near_ligand(pdb_path: str, ligand_atom_positions: list, cutoff: float = 5.0) -> set:
    """Which protein chains have at least one atom within cutoff of the
    kept ligand instance. Restricting to the ligand's own nominal chain
    (as an earlier version of this function did) is wrong for multi-chain
    functional units -- e.g. HIV-1 protease is an obligate homodimer
    where both chains form the single active site the ligand sits in;
    keeping only the ligand's own chain deletes half the binding site and
    leaves PDBFixer/OpenMM unable to build a sane Topology from what's
    left. Distance-based chain selection handles both that case and the
    original motivating one (a structure with multiple complex copies in
    the asymmetric unit, where only the copy near this ligand should be
    kept)."""
    import numpy as np
    lig_pos = np.array(ligand_atom_positions)

    chains = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if line[17:20].strip() in ('HOH',):
                continue
            try:
                pos = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if np.min(np.linalg.norm(lig_pos - pos, axis=1)) <= cutoff:
                chains.add(line[21])
    return chains


def prep_structure(pdb_id: str, pdb_path: str, ligand_resname: str, out_dir: str,
                    uniprot_acc: str = None) -> dict:
    ligand_raw = f"{out_dir}/{pdb_id}_ligand_raw.pdb"
    ligand_sdf = f"{out_dir}/{pdb_id}_ligand_h.sdf"
    receptor_input = f"{out_dir}/{pdb_id}_receptor_input.pdb"
    receptor_pkl = f"{out_dir}/{pdb_id}_receptor_h.pkl"

    n_atoms, chain_id = extract_single_ligand_instance(pdb_path, ligand_resname, ligand_raw)
    if n_atoms == 0:
        raise ValueError(f"No HETATM records found for resname {ligand_resname} in {pdb_path}")

    subprocess.run(f"obabel {ligand_raw} -O {ligand_sdf} -h", shell=True, check=True, capture_output=True)

    ligand_positions = []
    with open(ligand_raw) as f:
        for line in f:
            ligand_positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    keep_chains = chains_near_ligand(pdb_path, ligand_positions)

    # Restrict the receptor to chains actually near this ligand instance
    # (protein chains forming its binding site, not solvent/other-copy
    # chains far away), rather than the ligand's own nominal chain alone.
    with open(pdb_path) as fin, open(receptor_input, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                if line[21] in keep_chains:
                    fout.write(line)
            else:
                fout.write(line)

    if uniprot_acc:
        # Author residue numbers are not comparable across independent
        # PDB depositions (different construct boundaries/isoforms can
        # give the same author number to different physical residues --
        # see posegate.residue_mapping's docstring for the ERalpha case
        # that motivated this). Remap onto UniProt numbering, which is
        # consistent across every structure of the same protein, before
        # this receptor is ever mined.
        residue_map = build_residue_map(pdb_id, uniprot_acc)
        remapped_path = f"{out_dir}/{pdb_id}_receptor_input_remapped.pdb"
        n_remapped = remap_residue_numbers(receptor_input, remapped_path, residue_map)
        if n_remapped == 0:
            raise ValueError(f"{pdb_id}: SIFTS mapping to {uniprot_acc} matched 0 receptor atoms "
                              f"-- chain ID mismatch or wrong accession, refusing to mine unmapped")
        receptor_input = remapped_path

    prepare_receptor_pickle(receptor_input, receptor_pkl)

    return {'pdb_id': pdb_id, 'ligand_sdf': ligand_sdf, 'receptor_pdb': receptor_pkl}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, help="JSON file: [{pdb_id, pdb_path, ligand_resname}, ...]")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--out_manifest", required=True, help="Where to write the prepped-structure manifest")
    args = parser.parse_args()

    with open(args.manifest) as f:
        entries = json.load(f)

    prepped = []
    for entry in entries:
        print(f"Prepping {entry['pdb_id']} (ligand {entry['ligand_resname']})...")
        try:
            prepped.append(prep_structure(entry['pdb_id'], entry['pdb_path'], entry['ligand_resname'], args.out_dir))
        except Exception as e:
            print(f"  FAILED: {entry['pdb_id']} ({e})")

    with open(args.out_manifest, 'w') as f:
        json.dump(prepped, f, indent=2)
    print(f"Prepped {len(prepped)}/{len(entries)} structures. Manifest: {args.out_manifest}")


if __name__ == "__main__":
    main()
