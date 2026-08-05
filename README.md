<div align="center">

# posegate

**Mine the conserved contacts of a target from its own PDB structures, then check docked poses against them.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

</div>

Interaction-aware pose analysis needs a target-specific pharmacophore to check against, and for a
less-studied target no one has written that pharmacophore down. `posegate`'s main component
derives it from the target's own co-crystal structures: point it at several PDB entries for one
target, each with a different bound ligand, and it reports which receptor contacts recur and how
often. A second component applies that constraint to docked poses, reporting the interactions
each pose forms and whether it satisfies the mined contact. Interaction detection uses
[ProLIF](https://github.com/chemosim-lab/ProLIF); docking uses
[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

The miner is the validated part of this tool (five protein families, cross-checked against PLIP
and against visGReMLIN's published benchmark). The pose-ranking score is exploratory: it is a
modest improvement on raw docking score for the target it was fitted on, and its weights do not
transfer between targets. See [Validation](#validation) for both results in full.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Validation](#validation)
- [Design notes](#design-notes)
- [Related work](#related-work)
- [Project structure](#project-structure)
- [Testing](#testing)
- [Citation](#citation)
- [License](#license)

## Features

- **Conserved-contact mining** (`posegate.conserved_contacts`) — the main component. Given several
  co-crystal structures of one target, each bound to a different ligand, reports which receptor
  contacts recur across ligands and at what frequency. This yields a pharmacophore constraint for
  targets where the literature does not already provide one.
- **Pose autopsy** (`posegate.autopsy`) — steric clashes, hydrogen-bond geometry, aromatic
  contacts, and a checkable pharmacophore constraint (a named conserved H-bond, typically the one
  the miner found), combined into a fitted score for ranking candidates within a batch. The
  per-pose interaction report is reliable; the fitted ranking score is exploratory (see
  [Validation](#validation)).
- **Docking orchestration** (`posegate.docking`) — AutoDock Vina wrapper with ligand-size-aware
  search boxes and restraint-guided pose selection (pick the best-scoring pose among several
  that actually satisfies a required contact, not just the single top-ranked one).
- **Receptor preparation** (`posegate.receptor_prep`) — builds receptor bonds directly from
  PDBFixer/OpenMM's own `Topology`, rather than re-guessing them from PDB text (see
  [Design notes](#receptor-preparation) below for why that matters).

## Installation

```bash
conda env create -f environment.yml
conda activate posegate
pip install -e .
```

## Quickstart

### Mine conserved contacts across a PDB ensemble

Start here: the residue this reports is the pharmacophore constraint the other steps check
against.

```bash
# 1. Prepare each structure (ligand SDF + receptor pickle) from a manifest of
#    {pdb_id, pdb_path, ligand_resname} entries:
python scripts/prep_ensemble.py \
    --manifest my_targets.json --out_dir data/my_target --out_manifest data/my_target/prepped.json

# 2. Mine which receptor contacts are conserved across the ensemble:
python scripts/run_conserved_contact_miner.py --manifest data/my_target/prepped.json
```

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

### Batch dock and triage a compound library

```bash
python scripts/batch_dock.py \
    --receptor_pdb receptor.pdb --ligands_csv compounds.csv \
    --center 12.3 45.6 7.8 --exhaustiveness 32
```

## Validation

### Conserved-contact miner (five targets)

The miner was checked against independently verifiable ground truth and cross-validated against
[PLIP](https://doi.org/10.1093/nar/gkv315), a separately implemented interaction detector, on
five protein families spanning unrelated folds and binding chemistries.

| Family | Fold | Literature pharmacophore | Recovered by posegate | Agreed by PLIP |
|---|---|---|:---:|:---:|
| BRD4 | Bromodomain reader | Asn140 acetyl-lysine mimetic H-bond | ✅ | ✅ |
| CDK2 | Kinase | Leu83 hinge H-bond | ✅ | ✅ |
| Estrogen receptor alpha | Nuclear hormone receptor | Glu353/Arg394 "charge clamp" | ✅ | ✅ |
| HIV-1 protease | Homodimeric aspartic protease | Asp25/Asp25' catalytic dyad (both chains) | ✅ | ✅ |
| Carbonic anhydrase | Zinc metalloenzyme | Thr199 gatekeeper H-bond | ✅ | ✅ |

The literature-validated pharmacophore residue was recovered, and agreed on by both tools, in all
five families. Agreement across all conserved residues (≥50% ensemble frequency) is recall 0.32
and precision 0.71 against PLIP's broader output, with disagreement concentrated in PLIP's longer
tail of looser hydrophobic and water-bridge calls rather than in the pharmacophore residues.
Reproduce with `scripts/plip_ensemble_miner.py` and `scripts/compare_miners.py`.

We traced two representative disagreements to their causes, and neither is a detection error:

- **A CDK2 lysine contact that `posegate` reports and PLIP does not.** The donor-H···acceptor
  angle is 126.8°, just below ProLIF's default 130–180° acceptance window, on an otherwise
  genuine close contact. The two tools draw the hydrogen-bond angle cutoff in slightly different
  places, and this contact falls between them.
- **An apparent chain asymmetry in the HIV-1 protease catalytic aspartates.** Asp25 and Asp25'
  are not detected symmetrically in every structure. This reflects a real per-structure asymmetry
  in the deposited coordinates rather than a bug: pseudo-symmetric inhibitors are known to bind
  this homodimer with one aspartate engaged more closely than the other. Aggregated across the
  ensemble, the miner recovers the dyad from both chains.

### Pose-ranking score (five targets, exploratory)

This component is reported as a characterization of the approach's limits, not as a validated
screening method. Each target has its own DUD-E-style property-matched benchmark and its own
weights, fitted independently by `scripts/recalibrate_weights.py` (L2-regularized logistic
regression, not hand-tuned) against that target's own mined pharmacophore:

| Target | N | actives | raw Vina AUC | fitted AUC |
|---|---|---|---|---|
| CDK2 | 64 | 22 | 0.52 | 0.52 |
| BRD4 | 90 | 30 | 0.60 | 0.62 |
| Estrogen receptor alpha | 223 | 75 | 0.72 | 0.85 |
| HIV-1 protease | 225 | 75 | 0.72 | 0.77 |
| Carbonic anhydrase | 261 | 87 | 0.25 | 0.79 |

Fitted weight *magnitudes* are not comparable across targets at these sample sizes: with 64 to
261 compounds per target the confidence intervals overlap too much for a pairwise AUC comparison
to mean anything. What is comparable is each feature's *direction*. `scripts/compare_feature_weights.py`
fits each target's weights on standardized features and bootstraps 200 resamples per target,
reporting a feature's direction only when its sign is stable in at least 90% of resamples:

```
feature                       cdk2            brd4         eralpha             hiv              ca
----------------------------------------------------------------------------------------------------
conserved_hbond            +0.43           +0.73           +0.98           +0.45           +0.46
hbond_count                 unstable        -0.54           +0.30           -0.35            unstable
vina_score                   unstable        -0.61           -1.11           -0.49           +0.82
```

`conserved_hbond` is the one feature whose direction holds across all five targets (bootstrap
sign stability 0.94–1.00). Generic H-bond count reverses between targets, as does the raw Vina
score itself. This is consistent with the mined conserved contact, rather than the fitted score,
being the component of `posegate` that generalizes across targets — the closest this project
comes to a supported answer on that question, run on real docking rather than assumed.

**Carbonic anhydrase is worth a closer look, because its raw-Vina baseline of 0.25 means Vina**
**ranked actives *worse* than decoys on average, yet its mined Thr199 contact still discriminated**
**correctly.** We checked whether this was a bug in our own zinc handling before writing it up.
Docked poses place carbonic anhydrase inhibitors 1.4–2.9 Å from the catalytic zinc, consistent
with genuine coordination, against 1.8–4.2 Å for decoys — so pose selection is finding the
correct site. The zinc carries zero partial charge in the prepared receptor, and
[AutoDock Vina's own documentation](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html)
states that it disregards atomic charges on metal ions during scoring; the
[AutoDock4Zn](https://doi.org/10.1021/ci500209e) extension exists specifically to address this,
and we are not using it. Vina's raw score is unreliable here by a documented limitation of the
scoring function, not a defect in `posegate.receptor_prep`, and it is a clean illustration of why
a feature that doesn't depend on that score is useful.

Two limitations bound how far these results should be read. First, property-matched decoys carry
[analogue and decoy bias that a fitted model can learn in place of genuine protein-ligand
interaction](https://doi.org/10.1371/journal.pone.0220113); the fitted AUCs above are an upper
bound until validated against a bias-corrected benchmark such as
[LIT-PCBA](https://doi.org/10.1021/acs.jcim.0c00155). This bounds the fitted weights only, not
the miner, which uses no decoys. Second, per-target variance in raw docking performance (0.25 to
0.72 here) is expected rather than a sign of a broken setup — on LIT-PCBA even consensus scoring
fails to enrich on several of fifteen targets. Full discussion, including the metal-retention and
class-imbalance fixes this required, in `paper.md`.

## Design notes

### Docking and ranking

- **Ligand-size-aware search box.** The box is sized to each ligand's own conformer extent plus a
  fixed margin. A box substantially smaller than the ligand produces severe clashes, so a single
  fixed box size across a chemically diverse library is unsafe.
- **Restraint-guided pose selection.** Vina's scoring function has no restraint term, so it
  cannot be told to prefer poses making a particular contact. Rather than accepting Vina's single
  top-ranked pose, `posegate` requests several candidates and selects the best-scoring one that
  satisfies the required contact, falling back to the top-ranked pose when none does.
- **Percentile ranking within a batch.** A fixed absolute score threshold is only meaningful
  relative to the score distribution a given receptor and setup actually produces, so batch
  screening ranks candidates by percentile within the screened set instead.

### Receptor preparation

RDKit's native PDB bond perception is unreliable for full multi-residue proteins. With its
default proximity-bonding heuristic it can introduce spurious bonds at tight turns, and with that
heuristic disabled its residue-template matcher can leave most of a chain unbonded without
reporting an error. `posegate.receptor_prep` avoids both failure modes by building the receptor's
RDKit `Mol` from PDBFixer/OpenMM's `Topology.bonds()`, which OpenMM computes in order to run MD
simulations, rather than inferring connectivity from PDB text.

## Related work

`posegate` occupies a narrow scope relative to a few related tools:

- **[PoseBusters](https://github.com/maabuu/posebusters)** checks a pose's chemical and physical
  plausibility (bond lengths, ring planarity, stereochemistry, clashes) but not target-specific
  interaction recovery. The two are complementary: PoseBusters asks whether a pose is chemically
  valid, `posegate` asks which interactions a valid pose forms with the target.
- **[Errington et al. 2024/2025](https://doi.org/10.1186/s13321-025-01011-6)** measures how
  closely a predicted pose reproduces the interactions of a known reference crystal pose, also
  via ProLIF. That requires ground truth for comparison, whereas `posegate` addresses the more
  common screening case in which no reference pose exists.
- **[ParaDockS](https://github.com/cbaldauf/paradocks)** proposed essentially the same
  conserved-contact-mining idea over a decade ago, but its repository has had no commits since
  2015 and is not a maintained, installable package. `posegate` implements the idea on a
  currently maintained interaction-fingerprinting library (ProLIF).
- **[visGReMLIN](https://doi.org/10.1186/s12859-020-3347-7)** is the closest direct comparator to
  the conserved-contact miner. It takes the same input (an ensemble of one target's structures,
  each with a different ligand) but uses a different method: graph-mining of conserved 3D motifs
  rather than per-residue frequency aggregation over ProLIF fingerprints. It was released as a
  web server
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
  | `posegate`, 22-structure, specific (non-VdW) contacts only | 8/9 | 1/9 |

  The 9/9 indicates that the residue-level method recovers the same published binding site at
  lower cost. It is not evidence of better resolution: visGReMLIN's score is atom-level over 73
  complexes, and `posegate` produces no atom-level output to score on that denominator. On
  specific interactions both methods reach 8/9 and both miss the same residue, HIS84, which
  `posegate` recovers only as van der Waals proximity. PHE82, which visGReMLIN recovers through
  an aromatic motif, is recovered here as a hydrophobic contact rather than as pi-stacking.

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
  the table. Structure preparation is not bit-for-bit reproducible, since hydrogen placement
  varies slightly between runs, so individual `n_structures` counts can shift by one. The
  reference-residue scores above have been stable across reruns.
- **[FTMap](https://ftmap.bu.edu/)** / **[Fragment Hotspot Maps](https://fragment-hotspot-maps.ccdc.cam.ac.uk/)**
  identify druggable hot spots using FFT-accelerated probe-docking simulation and a statistical
  model built from the licensed Cambridge Structural Database respectively. Both need
  substantially heavier infrastructure and locate hot spots on a single structure. `posegate`'s
  miner works from structures already deposited in the public PDB and asks a narrower question:
  across several ligands known to bind this target, which receptor contacts recur?

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
├── data/                      # ensemble manifests; fetched structures are gitignored
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
