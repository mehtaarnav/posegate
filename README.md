# posegate

A computational chemistry toolkit for triaging and analyzing molecular docking poses.

## Overview

`posegate` provides utilities for:

- **Docking** (`posegate.docking`) — running and orchestrating docking workflows.
- **Ensemble analysis** (`posegate.ensemble`) — working with ensembles of docked poses.
- **Pose autopsy** (`posegate.autopsy`) — diagnosing docking poses via ProLIF-based steric clash
  detection, hydrogen bond geometry, aromatic contacts, and report generation.
- **Conserved-contact mining** (`posegate.conserved_contacts`) — mines a PDB ensemble of
  structures for one target (each bound to a different ligand) to automatically surface which
  receptor contacts are conserved across ligands, rather than requiring a hand-picked literature
  pharmacophore.
- **Utilities** (`posegate.utils`) — shared geometry and I/O helpers.

## Related work

`posegate`'s scope is deliberately narrow relative to a few related tools:

- **[PoseBusters](https://github.com/maabuu/posebusters)** checks a pose's chemical/physical
  plausibility (bond lengths, ring planarity, stereochemistry, clashes) but not target-specific
  interaction recovery. Complementary, not overlapping: PoseBusters asks "is this pose chemically
  valid?", `posegate` asks "does this valid pose engage the right interactions for this target?"
- **[Errington et al. 2024/2025](https://doi.org/10.1186/s13321-025-01011-6)** measures how well a
  predicted pose reproduces a *known reference crystal pose*'s interactions (also via ProLIF). That
  requires a ground-truth answer to compare against. `posegate` targets the opposite, more common
  screening situation: triaging candidates for which no reference pose exists.
- **[ParaDockS](https://github.com/cbaldauf/paradocks)** proposed essentially the same
  conserved-contact-mining idea over a decade ago, but its repository has had no commits since
  2015 — it's not a maintained, installable package. `posegate` re-implements the idea on a
  currently maintained interaction-fingerprinting library (ProLIF).
- **[FTMap](https://ftmap.bu.edu/)** / **[Fragment Hotspot Maps](https://fragment-hotspot-maps.ccdc.cam.ac.uk/)**
  identify druggable hot spots via, respectively, FFT-accelerated probe-docking simulation and a
  statistical model built from the licensed Cambridge Structural Database — both substantially
  heavier infrastructure aimed at finding hot spots on one structure. `posegate`'s miner instead
  asks a lighter-weight question of structures that already exist in the public PDB: across
  several different ligands already known to bind this target, which receptor contacts recur?
- **[visGReMLIN](https://doi.org/10.1186/s12859-020-3347-7)** is the closest direct comparator to
  the conserved-contact miner — same input specification (an ensemble of one target's structures,
  each with a different ligand), different method (graph-mining of conserved 3D motifs vs.
  per-residue frequency aggregation over ProLIF fingerprints). We have not yet run a head-to-head
  comparison against it specifically; see `paper.md`.

## Validation status

- **Interaction detection**: validated against real, independently-verifiable ground truth, and
  cross-checked against [PLIP](https://doi.org/10.1093/nar/gkv315) (an independent,
  differently-implemented interaction detector) on 5 protein families spanning unrelated folds —
  BRD4 (bromodomain), CDK2 (kinase), estrogen receptor alpha (nuclear receptor), HIV-1 protease
  (homodimeric aspartic protease), and carbonic anhydrase (zinc metalloenzyme), using 4-6
  co-crystal structures per family (`scripts/plip_ensemble_miner.py`, `scripts/compare_miners.py`).
  **The literature-validated pharmacophore residue was recovered, and agreed on by both tools, in
  all 5 families**: BRD4 Asn140, CDK2 Leu83 hinge H-bond, ERα Glu353/Arg394 charge clamp, HIV-1
  protease Asp25/Asp25' catalytic dyad (both chains), carbonic anhydrase Thr199 gatekeeper.
  Aggregate agreement across all conserved residues (≥50% ensemble frequency) is recall 0.32,
  precision 0.71 against PLIP's broader output — disagreement is concentrated in PLIP's longer
  tail of looser hydrophobic/water-bridge calls, not in the core pharmacophore signal. Full
  per-family numbers in `paper.md`.
- **Screening/ranking performance**: on a 65-compound BRD4 benchmark (22 actives, 43
  property-matched decoys), raw Vina score achieves AUC-ROC 0.53 (95% CI [0.37, 0.68]). The
  `posegate_score` weights (see `scripts/recalibrate_weights.py`) are fit by L2-regularized
  logistic regression against this benchmark's labels rather than hand-picked; honest
  cross-validated AUC-ROC is 0.62 — a real improvement, but calibrated on one target/65 compounds,
  not yet validated as general-purpose across targets. See `paper.md` for the full picture,
  including why generic H-bond count ends up *penalized* (property-matched decoys form just as
  many incidental H-bonds; only the specific conserved contact discriminates).

## Installation

See [Setup](#setup) below.

## Setup

```bash
conda env create -f environment.yml
conda activate posegate
pip install -e .
```

## Testing

```bash
pytest
```

## Project Structure

```
posegate/
├── posegate/       # library source
├── tests/          # unit tests
├── data/           # local data (gitignored)
└── scripts/        # entry-point scripts
```

## License

MIT
