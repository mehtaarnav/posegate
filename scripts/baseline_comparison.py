# posegate/scripts/baseline_comparison.py
"""Does the pipeline actually beat something dumber?

Everything this project has validated so far establishes that the
conserved-contact miner RECOVERS known selectivity-determining residues
across five protein families. None of it establishes that the machinery
doing the recovering is necessary: the literature-confirmed residues are,
by definition, binding-site residues, and a naive "list the residues near
the ligand" would return binding-site residues too. If a dumb baseline
recovers the same targets equally well, most of this apparatus is
decoration and the honest scope of the project shrinks accordingly.

Three baselines of increasing sophistication, each isolating one layer
the full pipeline adds:

  B1  geometry, single structure   -- residues with any atom within
      GEOM_CUTOFF of any ligand atom, ranked by number of contacting
      atom pairs. No interaction typing, no ensemble. The dumbest thing
      that could work.
  B2  ProLIF, single structure     -- that one structure's specific
      (non-VdW) interaction residues, ranked by how many distinct
      interaction types each forms. Adds interaction typing, still no
      ensemble.
  B3  geometry, ensemble           -- residues ranked by the FRACTION of
      structures in which they fall within GEOM_CUTOFF. Adds the
      ensemble/conservation logic, still no interaction typing.

  FULL  the actual pipeline (mine_conserved_contacts): ProLIF interaction
      typing + ensemble conservation frequency, read from each target's
      already-computed mined_result.json.

B1 and B2 are single-structure methods, so running them on one
arbitrarily-chosen member would make the result depend on that choice.
Each is instead run on EVERY structure in the ensemble and reported as a
mean recovery rate across those runs -- which also directly measures how
much the ensemble is worth, since B3/FULL get to see all of them.

Metric is recovery@10: does the literature-confirmed selectivity residue
for that protein appear in the method's top 10 residues? Those targets
(TARGETS below) are the ones independently confirmed against published
literature earlier in this project -- not chosen here, and not chosen
with any knowledge of how the baselines would score.
"""

import glob
import json
import os
import sys

import numpy as np
from rdkit import Chem

sys.path.insert(0, os.path.dirname(__file__))
from posegate.autopsy import build_ifp
from posegate.conserved_contacts import SPECIFIC_INTERACTIONS, _load_structure

GEOM_CUTOFF = 4.5
TOP_K = 10

# Literature-confirmed selectivity residues, in each protein's own mined
# (SIFTS-remapped UniProt) numbering. Every one was verified against
# published sources earlier in this project -- see the HOLDOUT_RESULT_*
# files and ca/cdk_family_selectivity.py docstrings. A protein counts as
# "recovered" if ANY of its listed targets makes the method's top-10.
#
# 'preregistered' marks HOW the target was chosen, which matters more
# than it might appear. Four of these (CA II, CA IX, CDK2, CDK9) were
# identified by reading the pipeline's own top-10 output and THEN
# verified against literature. For those, FULL scoring a hit is
# circular -- the residue is in FULL's top-10 by construction, because
# that is where it was found. Only the four marked preregistered=True
# were predicted from published literature BEFORE any mining ran on that
# family (see PREREGISTRATION_*.md), so only those measure anything
# about FULL. Both subsets are reported separately below; the
# preregistered subset is the one that carries evidential weight.
TARGETS = {
    'ca_verified':      {'label': 'CA II (P00918)',    'targets': [130],      'preregistered': False},
    'ca9_verified':     {'label': 'CA IX (Q16790)',    'targets': [262, 263], 'preregistered': False},
    'cdk2_verified':    {'label': 'CDK2 (P24941)',     'targets': [83],       'preregistered': False},
    'cdk9_verified':    {'label': 'CDK9 (P50750)',     'targets': [106],      'preregistered': False},
    'trypsin_verified2': {'label': 'Trypsin (P00760)',  'targets': [194],      'preregistered': True},
    'ache_verified':    {'label': 'AChE (P22303)',     'targets': [326, 328], 'preregistered': True},
    'cox1_verified':    {'label': 'COX-1 (P05979)',    'targets': [523],      'preregistered': True},
    'cox2_verified':    {'label': 'COX-2 (P35354)',    'targets': [509],      'preregistered': True},
}


def load_ensemble(data_dir):
    """Reconstructs the prepped-structure list from an output directory,
    pairing each receptor pickle with its ligand SDF by PDB ID."""
    structures = []
    for pkl in sorted(glob.glob(os.path.join(data_dir, '*_receptor_h.pkl'))):
        pdb_id = os.path.basename(pkl).replace('_receptor_h.pkl', '')
        sdf = os.path.join(data_dir, f'{pdb_id}_ligand_h.sdf')
        if os.path.exists(sdf):
            structures.append({'pdb_id': pdb_id, 'ligand_sdf': sdf, 'receptor_pdb': pkl})
    return structures


def _residue_number(atom):
    info = atom.GetPDBResidueInfo()
    return info.GetResidueNumber() if info is not None else None


def geometric_contacts(structure, cutoff=GEOM_CUTOFF):
    """{residue_number: n_contacting_atom_pairs} by raw distance, no
    interaction typing. None if the structure fails to load."""
    lig, rec = _load_structure(structure)
    if lig is None or rec is None:
        return None
    lig_pos = np.array([lig.GetConformer().GetAtomPosition(i) for i in range(lig.GetNumAtoms())])
    rec_conf = rec.GetConformer()

    counts = {}
    for atom in rec.GetAtoms():
        resnum = _residue_number(atom)
        if resnum is None:
            continue
        pos = np.array(rec_conf.GetAtomPosition(atom.GetIdx()))
        n_close = int(np.sum(np.linalg.norm(lig_pos - pos, axis=1) <= cutoff))
        if n_close:
            counts[resnum] = counts.get(resnum, 0) + n_close
    return counts


def prolif_contacts(structure):
    """{residue_number: n_distinct_specific_interaction_types} for one
    structure. None if it fails to load or ProLIF errors on it."""
    lig, rec = _load_structure(structure)
    if lig is None or rec is None:
        return None
    try:
        ifp = build_ifp(lig, rec)
    except Exception:
        return None
    counts = {}
    for (_, pres), interactions in ifp.items():
        specific = SPECIFIC_INTERACTIONS.intersection(interactions)
        if not specific:
            continue
        digits = ''.join(c for c in str(pres) if c.isdigit())
        if digits:
            counts[int(digits)] = len(specific)
    return counts


def top_k(counts, k=TOP_K):
    return [r for r, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:k]]


def full_pipeline_top_k(data_dir, k=TOP_K):
    """Top-k distinct residues from the already-computed mined result,
    VdWContact excluded -- the same convention the rest of the project
    uses for reading its own output."""
    with open(os.path.join(data_dir, 'mined_result.json')) as f:
        data = json.load(f)
    seen = []
    for row in data['mined']:
        if row['interaction'] == 'VdWContact':
            continue
        digits = ''.join(c for c in row['residue'] if c.isdigit())
        if not digits:
            continue
        n = int(digits)
        if n not in seen:
            seen.append(n)
        if len(seen) >= k:
            break
    return seen


def main():
    print(f"Baseline comparison: recovery@{TOP_K} of literature-confirmed "
          f"selectivity residues\n")
    header = (f"{'Target':<20}{'B1 geom/1':>11}{'B2 prolif/1':>13}"
              f"{'B3 geom/ens':>13}{'FULL':>8}   targets")
    print(header)
    print('-' * len(header))

    totals = {'B1': [], 'B2': [], 'B3': [], 'FULL': []}
    prereg = {'B1': [], 'B2': [], 'B3': [], 'FULL': []}

    for data_dir, meta in TARGETS.items():
        path = os.path.join('data', data_dir)
        structures = load_ensemble(path)
        targets = set(meta['targets'])
        if not structures:
            print(f"{meta['label']:<20}  (no prepped structures found)")
            continue

        # B1/B2: single-structure methods, run on every member so the
        # score doesn't depend on an arbitrary choice of which one.
        b1_hits, b1_n = 0, 0
        b2_hits, b2_n = 0, 0
        per_structure_geom = []
        for s in structures:
            geom = geometric_contacts(s)
            if geom is not None:
                per_structure_geom.append(set(geom.keys()))
                b1_n += 1
                if targets & set(top_k(geom)):
                    b1_hits += 1
            prolif = prolif_contacts(s)
            if prolif is not None:
                b2_n += 1
                if targets & set(top_k(prolif)):
                    b2_hits += 1

        # B3: ensemble geometry -- fraction of structures in which each
        # residue falls inside the cutoff.
        ens_counts = {}
        for residues in per_structure_geom:
            for r in residues:
                ens_counts[r] = ens_counts.get(r, 0) + 1
        b3_hit = bool(targets & set(top_k(ens_counts))) if ens_counts else False

        full_hit = bool(targets & set(full_pipeline_top_k(path)))

        b1_rate = b1_hits / b1_n if b1_n else 0.0
        b2_rate = b2_hits / b2_n if b2_n else 0.0
        for bucket, val in (('B1', b1_rate), ('B2', b2_rate),
                             ('B3', 1.0 if b3_hit else 0.0),
                             ('FULL', 1.0 if full_hit else 0.0)):
            totals[bucket].append(val)
            if meta['preregistered']:
                prereg[bucket].append(val)

        flag = ' [prereg]' if meta['preregistered'] else ''
        print(f"{meta['label']:<20}{b1_rate:>10.0%}{b2_rate:>13.0%}"
              f"{'YES' if b3_hit else 'no':>13}{'YES' if full_hit else 'no':>8}   "
              f"{sorted(targets)}{flag}")

    print('-' * len(header))
    print(f"{'MEAN (all 8)':<20}{np.mean(totals['B1']):>10.0%}{np.mean(totals['B2']):>13.0%}"
          f"{np.mean(totals['B3']):>12.0%}{np.mean(totals['FULL']):>8.0%}")
    print(f"{'MEAN (prereg only)':<20}{np.mean(prereg['B1']):>10.0%}{np.mean(prereg['B2']):>13.0%}"
          f"{np.mean(prereg['B3']):>12.0%}{np.mean(prereg['FULL']):>8.0%}")
    print("\nB1/B2 are per-structure rates averaged over every ensemble member; "
          "B3/FULL see the whole ensemble at once, so they are hit/miss per target.")
    print("The 'all 8' FULL column is inflated: 4 of those targets were found BY "
          "reading FULL's own output, so it cannot miss them. Only the 4 marked "
          "[prereg] were predicted from literature before mining ran -- that row "
          "is the one that measures anything.")


if __name__ == "__main__":
    main()
