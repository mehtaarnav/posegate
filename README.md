<div align="center">

# posegate

**Interaction-aware triage for docked protein-ligand poses, and automated conserved-contact mining across PDB ensembles.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

</div>

Docking gives you a pose and a score. `posegate` tells you *why* that pose is or isn't
trustworthy — which interactions it makes, whether they include the residues a real
inhibitor is known to engage, and how consistent that pharmacophore is across the target's
whole structural history — using [ProLIF](https://github.com/chemosim-lab/ProLIF) for
interaction detection and [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) for
docking.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Validation](#validation)
- [Related work](#related-work)
- [Project structure](#project-structure)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)

## Features

- **Pose autopsy** (`posegate.autopsy`) — steric clashes, hydrogen-bond geometry, aromatic
  contacts, and a checkable pharmacophore constraint (a named conserved H-bond), combined into
  a fitted score for ranking candidates within a batch.
- **Conserved-contact mining** (`posegate.conserved_contacts`) — given several real co-crystal
  structures of one target, each bound to a different ligand, automatically surfaces which
  receptor contacts recur across ligands. This is the data-driven equivalent of a hand-picked
  literature pharmacophore, without requiring you to already know the literature.
- **Docking orchestration** (`posegate.docking`) — AutoDock Vina wrapper with ligand-size-aware
  search boxes and restraint-guided pose selection (pick the best-scoring pose among several
  that actually satisfies a required contact, not just the single top-ranked one).
- **Receptor preparation** (`posegate.receptor_prep`) — builds receptor bonds directly from
  PDBFixer/OpenMM's own `Topology`, rather than re-guessing them from PDB text (see
  [Notes](#notes-on-receptor-preparation) below for why that matters).

## Installation

```bash
conda env create -f environment.yml
conda activate posegate
pip install -e .
```

## Quickstart

### Autopsy a single docked pose

```python
from posegate.autopsy import generate_autopsy_report

report = generate_autopsy_report(
    ligand_sdf_path="ligand_docked.sdf",
    receptor_pdb_path="receptor.pkl",  # from posegate.receptor_prep, see below
    vina_score=-8.4,
)
print(report["decision"])       # PRIORITIZE / REVIEW / REJECT
print(report["hbonds"])         # detected hydrogen bonds
print(report["conserved_hbond"])  # hit on the named pharmacophore constraint, if any
```

Or from the command line:

```bash
python scripts/run_autopsy.py --receptor receptor_h.pdb --ligand ligand_docked.sdf --score -8.4
```

### Mine conserved contacts across a PDB ensemble

```bash
# 1. Prepare each structure (ligand SDF + receptor pickle) from a manifest of
#    {pdb_id, pdb_path, ligand_resname} entries:
python scripts/prep_ensemble.py \
    --manifest my_targets.json --out_dir data/my_target --out_manifest data/my_target/prepped.json

# 2. Mine which receptor contacts are conserved across the ensemble:
python scripts/run_conserved_contact_miner.py --manifest data/my_target/prepped.json
```

### Batch dock and triage a compound library

```bash
python scripts/batch_dock.py \
    --receptor_pdb receptor.pdb --ligands_csv compounds.csv \
    --center 12.3 45.6 7.8 --exhaustiveness 32
```

## Validation

Both components were checked against real, independently-verifiable ground truth, and
cross-validated against [PLIP](https://doi.org/10.1093/nar/gkv315) — an established,
independently-implemented interaction detector — on five protein families spanning unrelated
folds and binding chemistries.

| Family | Fold | Literature pharmacophore | Recovered by posegate | Agreed by PLIP |
|---|---|---|:---:|:---:|
| BRD4 | Bromodomain reader | Asn140 acetyl-lysine mimetic H-bond | ✅ | ✅ |
| CDK2 | Kinase | Leu83 hinge H-bond | ✅ | ✅ |
| Estrogen receptor alpha | Nuclear hormone receptor | Glu353/Arg394 "charge clamp" | ✅ | ✅ |
| HIV-1 protease | Homodimeric aspartic protease | Asp25/Asp25' catalytic dyad (both chains) | ✅ | ✅ |
| Carbonic anhydrase | Zinc metalloenzyme | Thr199 gatekeeper H-bond | ✅ | ✅ |

**The literature-validated pharmacophore residue was recovered, and agreed on by both tools, in
all 5 families.** Aggregate agreement across all conserved residues (≥50% ensemble frequency) is
recall 0.32 / precision 0.71 against PLIP's broader output; disagreement is concentrated in
PLIP's longer tail of looser hydrophobic/water-bridge calls, not in the core pharmacophore
signal. Two representative disagreements were diagnosed and traced to genuine geometric
threshold effects and real per-structure biological asymmetry, not detection bugs — see
[`paper.md`](paper.md) for the full breakdown and reproduction scripts
(`scripts/plip_ensemble_miner.py`, `scripts/compare_miners.py`).

**Screening/ranking performance** is a separate, weaker result, reported honestly rather than
polished: on a 65-compound BRD4 benchmark (22 actives, 43 property-matched decoys), raw Vina
score achieves AUC-ROC 0.53 (95% CI [0.37, 0.68]) — indistinguishable from random at this
sample size. `posegate`'s fitted `posegate_score` (weights from L2-regularized logistic
regression against this benchmark, not hand-picked — see `scripts/recalibrate_weights.py`)
reaches a cross-validated AUC-ROC of 0.62 *on the target it was fit on*. Notably, generic H-bond
count gets a *penalizing* weight in the fit — property-matched decoys form just as many
incidental H-bonds as actives, so only the *specific* conserved contact actually discriminates.

**Transferability was tested directly, and the fixed weights do not transfer.** Applying the
same BRD4-fitted weights unmodified to an equivalent 51-compound CDK2 benchmark (CDK2's own
Leu83 hinge contact substituted for BRD4's Asn140) gives AUC-ROC 0.35 — worse than raw Vina's
already-weak 0.37 on that benchmark. The reason is diagnosed, not just observed: on CDK2, actives
average *more* H-bonds than decoys (2.0 vs. 1.4) — the opposite pattern from BRD4 — so the
BRD4-derived H-bond penalty punishes exactly the signal that's real on CDK2. The
conserved-contact-hit and clash-count features keep the correct-direction relationship to
activity on both targets; generic H-bond count does not. Practical upshot: whether a feature
predicts activity is itself target-dependent, so per-target refitting (or restricting to
mechanistically target-general features like the conserved-contact hit) is necessary, not
optional. Full discussion in `paper.md`.

### Notes on receptor preparation

RDKit's native PDB bond-perception is unreliable for full multi-residue proteins: with its
default proximity-bonding heuristic it can invent spurious bonds on tight turns, and with that
heuristic disabled its residue-template matcher can silently leave the vast majority of a chain
completely unbonded. `posegate.receptor_prep` sidesteps both failure modes by building the
receptor's RDKit `Mol` directly from PDBFixer/OpenMM's own `Topology.bonds()` — bonds OpenMM
already computed correctly to run MD simulations — rather than re-guessing connectivity from
written-out PDB text.

## Related work

`posegate`'s scope is deliberately narrow relative to a few related tools:

- **[PoseBusters](https://github.com/maabuu/posebusters)** checks a pose's chemical/physical
  plausibility (bond lengths, ring planarity, stereochemistry, clashes) but not target-specific
  interaction recovery. Complementary, not overlapping: PoseBusters asks "is this pose chemically
  valid?"; `posegate` asks "does this valid pose engage the right interactions for this target?"
- **[Errington et al. 2024/2025](https://doi.org/10.1186/s13321-025-01011-6)** measures how well a
  predicted pose reproduces a *known reference crystal pose*'s interactions (also via ProLIF).
  That requires a ground-truth answer to compare against. `posegate` targets the opposite, more
  common screening situation: triaging candidates for which no reference pose exists.
- **[ParaDockS](https://github.com/cbaldauf/paradocks)** proposed essentially the same
  conserved-contact-mining idea over a decade ago, but its repository has had no commits since
  2015 — it's not a maintained, installable package. `posegate` re-implements the idea on a
  currently maintained interaction-fingerprinting library (ProLIF).
- **[visGReMLIN](https://doi.org/10.1186/s12859-020-3347-7)** is the closest direct comparator to
  the conserved-contact miner — same input specification (an ensemble of one target's structures,
  each with a different ligand), different method (graph-mining of conserved 3D motifs vs.
  per-residue frequency aggregation over ProLIF fingerprints). It was released as a web server
  only, with no source code or package, and both advertised URLs are currently unreachable, so it
  cannot be run on new data. `scripts/compare_visgremlin.py` instead compares against its
  published result on its own CDK case study, which scores recovery of the Schonbrunn et al. CDK
  binding site (26 atoms across 9 residues):

  | | reference residues at the interface | of those, at freq ≥0.33 |
  |---|---|---|
  | visGReMLIN (73 complexes, atom-level motifs) | 8/9 — 18/26 atoms, 69% | n/a |
  | `posegate`, 6-structure ensemble | 9/9 | 8/9 |
  | `posegate`, 19-structure ensemble | 8/9 | 4/9 |
  | `posegate`, 22-structure ensemble (union) | 9/9 | 4/9 |
  | `posegate`, 22-structure, specific (non-VdW) contacts only | 7/9 | 1/9 |

  The 9/9 indicates that the residue-level method recovers the same published binding site at
  lower cost. It is not evidence of better resolution: visGReMLIN's score is atom-level over 73
  complexes, and `posegate` produces no atom-level output to score on that denominator. PHE82 is
  a real miss. visGReMLIN identified its aromatic contacts, while ProLIF registers them only as
  van der Waals proximity.

  Ensemble composition affects the result more than ensemble size does. Growing the ensemble from
  6 to 19 structures degraded it: ASP145 dropped out of the output entirely, because 14 of the
  added structures belong to one fragment-screen deposition series (6Q3B–6Q4K) that binds the
  hinge without reaching the DFG region. Restricting that series to its drug-like members (≥10
  heavy atoms) did not restore ASP145, so the cause is chemotype homogeneity rather than ligand
  size. Which residues appear at all is reasonably stable (8 to 9 of 9 across ensembles), but the
  frequencies are not; they describe how much of a given ensemble's chemistry touches each
  residue, and only describe the pocket when the ensemble spans distinct chemical series. Curate
  these ensembles for scaffold diversity rather than for size. Reproduce with:

  ```bash
  bash scripts/fetch_cdk2_candidates.sh
  ```
  ```bash
  python scripts/prep_ensemble.py --manifest data/ensemble_cdk2/cdk2_manifest22.json --out_dir data/ensemble_cdk2 --out_manifest data/ensemble_cdk2/cdk2_prepped22.json
  ```
  ```bash
  python scripts/run_conserved_contact_miner.py --manifest data/ensemble_cdk2/cdk2_prepped22.json --out_json data/ensemble_cdk2/cdk2_freq22.json
  ```
  ```bash
  python scripts/compare_visgremlin.py --freq_json data/ensemble_cdk2/cdk2_freq22.json
  ```

  Substitute `cdk2_manifest.json` (6 structures) or `cdk2_manifest19.json` for the other rows of
  the table. Structure preparation is not bit-for-bit deterministic — hydrogen placement varies
  slightly between runs — so individual `n_structures` counts can shift by one; the
  reference-residue scores above have been stable across reruns.
- **[FTMap](https://ftmap.bu.edu/)** / **[Fragment Hotspot Maps](https://fragment-hotspot-maps.ccdc.cam.ac.uk/)**
  identify druggable hot spots via, respectively, FFT-accelerated probe-docking simulation and a
  statistical model built from the licensed Cambridge Structural Database — both substantially
  heavier infrastructure aimed at finding hot spots on one structure. `posegate`'s miner instead
  asks a lighter-weight question of structures that already exist in the public PDB: across
  several different ligands already known to bind this target, which receptor contacts recur?

## Project structure

```
posegate/
├── posegate/                  # library source
│   ├── autopsy.py             # per-pose interaction detection + scoring
│   ├── conserved_contacts.py  # PDB-ensemble conserved-contact miner
│   ├── docking.py             # AutoDock Vina orchestration
│   └── receptor_prep.py       # PDBFixer/OpenMM-based receptor preparation
├── scripts/                   # CLI entry points (see Quickstart)
├── tests/                     # unit tests
├── data/                      # local data (gitignored)
├── paper.md / paper.bib       # JOSS submission draft
└── environment.yml            # conda environment spec
```

## Testing

```bash
pytest
```

## Citation

A JOSS submission is in preparation (`paper.md`). Until it's published, please cite this
repository directly.

## License

[MIT](LICENSE)
