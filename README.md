<div align="center">

# posegate

**Mine the conserved contacts of a target from its own PDB structures, compare them across a protein family, and check docked poses against them.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

</div>

Interaction-aware pose analysis needs a target-specific pharmacophore to check against, and for a
less-studied target nobody has written that pharmacophore down. `posegate` derives it from the
target's own co-crystal structures: point it at several PDB entries for one target, each with a
different bound ligand, and it reports which receptor contacts recur and how often, with a
leave-one-out self-validation score that needs no external ground truth. A second layer compares
those contacts across the isoforms of a protein family to surface candidate selectivity handles.
Interaction detection uses [ProLIF](https://github.com/chemosim-lab/ProLIF); docking uses
[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina).

## Status — read this first

This project has been tested harder than it has been promoted. Both are recorded below.

**Supported by evidence:**

> posegate systematises and ranks what the structural record already encodes — recovering
> literature-confirmed selectivity residues across several protein families, and doing so better
> than sequence variability or naive geometry alone.

**Tested and _not_ supported:**

> posegate's output predicts experimentally achievable isoform selectivity.

That second claim was never stated outright, but it is what a reader naturally infers from
"selectivity-mapping tool." It was tested twice — at pair level against ~43,000 ChEMBL Ki
measurements, and ligand-conditioned on compounds with solved structures — and failed both times,
with pre-specified controls. See [Negative results](#negative-results). Development stopped on
that thesis by a decision rule fixed before the tests were run.

The pose-ranking score is exploratory: a modest improvement on raw docking score for the target it
was fitted on, with weights that do not transfer between targets.

## Contents

- [Features](#features)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Validation](#validation)
- [Negative results](#negative-results)
- [Known limitations](#known-limitations)
- [Correctness history](#correctness-history)
- [Design notes](#design-notes)
- [Project structure](#project-structure)
- [Testing](#testing)
- [License](#license)

## Features

- **Conserved-contact mining** (`posegate.conserved_contacts`) — the core. Given several co-crystal
  structures of one target, each bound to a different ligand, reports which receptor contacts recur
  and at what frequency. Every frequency carries a 95% Wilson confidence interval, because a raw
  frequency alone does not distinguish a 6-structure ensemble from a 22-structure one.
- **Leave-one-out self-validation** — holds out each structure in turn, mines the rest, and checks
  whether the held-out ligand actually contacts the predicted residues. Produces a confidence
  number for *this* ensemble without any external ground truth, which is the situation the tool
  exists for. Reported automatically on every run, alongside an ensemble-size reliability tier
  derived from a measured accuracy/variance curve.
- **Cross-isoform selectivity mapping** (`posegate.selectivity`) — compares mined contacts across
  N isoforms of a family via pairwise sequence alignment, separating the pan-family catalytic
  scaffold from positions where residue identity diverges. `scripts/ca_family_matrix.py` renders
  this as a residue-by-isoform matrix ranked by variability.
- **Cross-structure residue identity** (`posegate.residue_mapping`) — remaps every structure onto
  UniProt numbering via SIFTS and rejects any structure whose SIFTS mapping does not match the
  declared accession. Not optional; see [Correctness history](#correctness-history).
- **Multi-coordinate numbering display** (`posegate.numbering_display`) — shows each mined residue's
  author numbering alongside its UniProt number, and flags inconsistency across the ensemble
  rather than hiding it.
- **Pose autopsy** (`posegate.autopsy`) — steric clashes, hydrogen-bond geometry, aromatic contacts,
  geometric metal-coordination checks, and a checkable pharmacophore constraint, combined into a
  fitted score. The per-pose interaction report is reliable; the fitted ranking score is exploratory.
- **Docking orchestration** (`posegate.docking`) — AutoDock Vina wrapper with ligand-size-aware
  search boxes and restraint-guided pose selection.
- **Receptor preparation** (`posegate.receptor_prep`) — builds receptor bonds from PDBFixer/OpenMM's
  own `Topology` rather than re-guessing them from PDB text.

## Installation

```bash
conda env create -f environment.yml
conda activate posegate
pip install -e .
```

`obabel` (Open Babel CLI) is a runtime dependency with no pip wheel bundling the executable. If it
is missing you get an explicit error naming it and the install command, not a bare exit code.

## Quickstart

### Mine a target end to end

One command: fetch, detect ligands, prepare, mine, self-validate.

```bash
python scripts/mine_target.py \
    --pdb_ids 4EY5 4EY6 4EY7 4M0E 4M0F \
    --out_dir data/my_target \
    --uniprot_acc P22303
```

`--uniprot_acc` is **required**. A keyword-built PDB ID list can silently include a different
protein — an "ERα" list once pulled in ERβ and ERR-γ structures — and independent depositions of
the same protein use inconsistent author numbering. Both corrupt mining without raising an error.

### Compare a whole protein family

```bash
python scripts/ca_family_matrix.py
```

Produces a residue-by-isoform matrix: invariant positions (shared catalytic scaffold, not
selectivity-exploitable) versus variable positions (candidate selectivity handles), ranked by how
many distinct amino acids occur across the family.

### Autopsy a docked pose

```python
from posegate.autopsy import generate_autopsy_report

report = generate_autopsy_report(
    ligand_sdf_path="ligand_docked.sdf",
    receptor_pdb_path="receptor.pkl",
    vina_score=-8.4,
)
print(report["decision"])         # PRIORITIZE / REVIEW / REJECT
print(report["conserved_hbond"])  # hit on the mined pharmacophore constraint
```

## Validation

### Conserved-contact recovery

The miner recovers textbook active-site chemistry across structurally unrelated folds — zinc
hydrolase, kinase, serine protease, cholinesterase, heme peroxidase, bromodomain, aspartyl protease:

| target | LOO top-1 | n | recovers |
|---|---|---|---|
| BRD4 | 100% | 11 | Asn140 (the literature-hardcoded contact this tool was built to replace), Ile146 |
| COX-1 | 83% | 12 | Ile523 |
| Trypsin | 81% | 16 | Asp189 S1 pocket (classical numbering) |
| ERα | 75% | 12 | Leu387/Leu391, Phe404, Glu353 |
| CDK2 | 73% | 15 | Leu83 hinge |
| CA II | 67% | 12 | His94, Thr198/199 gatekeeper |
| HIV protease | 40% | 10 | Asp25/Gly27 catalytic dyad |
| AChE | 39% | 13 | Phe295/Phe297 acyl pocket |

Every row was produced by the current pipeline. An earlier revision of this table reported HIV
protease at 87%, taken from a run that predated mandatory SIFTS remapping and accession
verification, on an ID list screened for neither mixed isolates nor resistance mutants. Re-run on
14 title-verified wild-type structures with remapping, it scores 40%. The higher figure was
inflated by an unverified ensemble; the drop is a correction, not a regression. Why this target
scores low — large symmetric homodimer site, chemically very heterogeneous inhibitor set — is
unexamined.

### Selectivity residues, confirmed against literature

Seven human carbonic anhydrase isoforms compared in one pass. All eight invariant positions are
textbook CA catalytic machinery. Of the variable positions:

| position (classical) | difference | status |
|---|---|---|
| 91 | Ile (CA2) / Phe (CA1) | **confirmed** — literature calls it *the* highest-variability position and a named selectivity "hot-spot" |
| 131 / 132 / 135 | Phe/Gly/Val (CA2) vs Val/Asp/Leu (CA9) vs Ala/Ser/Ser (CA12) | **confirmed** — all nine cells match a published subpocket design rule exactly |
| 200 | Thr (CA2) / His (CA1) | **confirmed** — swapping it shifts CA II anion-inhibition Ki toward CA I |
| 19 | — | **false positive**, reported |

Position 91 is the strongest result: the matrix's own ranking metric independently reproduced the
literature's characterisation of it as the most variable position, ranking it first of ten
unprompted. Position 200 could not have been found by the earlier three-isoform comparison, because
CA1 was not in it — widening the family surfaced a determinant a narrower view structurally could
not reach.

Two families were **pre-registered** — predictions committed to git before any structure search
(`PREREGISTRATION_*.md`) — and both confirmed.

### Against baselines

Four methods on the pre-registered subset (n=4), isolating what each layer contributes:

| method | interaction typing | ensemble | recovery@10 |
|---|---|---|---|
| geometry, single structure | no | no | 65% |
| ProLIF, single structure | yes | no | 41% |
| geometry, ensemble | no | yes | 50% |
| **full pipeline** | yes | yes | **100%** |

The components are complementary, not redundant. A useful negative finding: ensemble geometry
(50%) is *worse* than single-structure geometry (65%) — ranking by fraction-of-structures-in-contact
favours the always-occupied pocket core and demotes rim residues, where selectivity determinants
often sit.

Caveat stated plainly: n=4, so 4/4 against a 65% baseline is p ≈ 0.18. Directionally consistent,
not significant. And the dumbest baseline already gets two-thirds of these.

### Against a sequence-only baseline

Because the selectivity matrix ranks by a *sequence* property, a plain alignment reproducing the
result would make the structural pipeline decoration. Tested: position 91 ranks 11/260 by
whole-protein sequence variability, and the confirmed Thr200 discriminator ranks **159/260**. The
most sequence-variable positions protein-wide are overwhelmingly surface loops. The structural
filter does the discriminative work. Threat retired.

## Negative results

Recorded as prominently as the positive ones.

**Pair-level test** (`scripts/selectivity_vs_experiment.py`) — ~43,000 ChEMBL Ki measurements, 21
isoform pairs, one pre-specified statistic. Does contact divergence predict experimental
selectivity? ρ = +0.244, **p = 0.287**. Whole-protein control ρ = +0.199; 11.3% of random
non-contact position sets match or beat it. Post-hoc alternative observables are worse, not better.

**Why it failed** — variance decomposition over 38,151 paired measurements: **74.5% of selectivity
variance is within pairs (ligand-driven)**, only 25.5% between pairs. Collapsing each pair to a
median discarded three quarters of the signal by construction.

**Ligand-conditioned test** (`scripts/ligand_conditioned_test.py`) — the corrected version, using
compounds whose structures are in the mined ensembles so the exact contacted positions are known:

| | ρ | p |
|---|---|---|
| main — divergent contacts | −0.171 | 0.193 |
| control — invariant contacts | +0.341 | **0.0076** |

The pre-specified control succeeded and the hypothesis failed, pointing the wrong way. Compounds
touching *zero* divergent positions were the most selective in the sample.

**The deeper obstacle:** compounds with crystal structures are systematically the non-selective
ones. Pan-CA sulfonamides crystallise readily and have decades of study; genuinely isoform-selective
inhibitors are recent and largely structurally unsolved. The data needed to test this properly is
largely absent from the PDB.

Full detail in `EXPERIMENTAL_VALIDATION_RESULT.md`, `LIGAND_CONDITIONED_RESULT.md`, and
`THREATS_TO_CONCLUSION.md`.

## Known limitations

- **Asymmetric multi-chain receptors are unsupported.** Chymotrypsin-class targets — several
  non-identical chains from one cleaved polypeptide — can have chain letters assigned inconsistently
  across depositions. Detected and reported as an explicit unsupported-target-class warning rather
  than silently producing a plausible-looking wrong answer.
- **Multi-residue and peptide ligands** are not handled by ligand auto-detection.
- **Reliability tiers were derived from one target** (CDK2) and applied to all.
- **LOO measures internal consistency, not correctness.** A homogeneous ensemble from a single
  medicinal-chemistry series scores high and means little.
- **Top-N cutoffs are arbitrary.** One genuine member of the CA selective pocket (position 67) was
  missed at top-10.
- **Scope is enzymes.** No GPCRs, ion channels, or protein–protein interfaces.

## Correctness history

Several bug classes were found by evidence rather than by inspection, each surfacing as an
implausible result that was then traced. Recorded because the failure modes generalise.

| symptom | cause | fix |
|---|---|---|
| ERα LOO top-1 **0%** | independent depositions assign the same author residue number to different physical residues; the "ERα" ensemble also contained ERβ and ERR-γ | mandatory SIFTS remapping + accession verification (**0% → 80%**) |
| CA13 LOO top-1 **0%** at HIGH reliability | residue labels carried chain letters, so equivalent crystallographic copies counted as separate residues (`PHE132.A` 0.33 + `PHE132.B` 0.40 instead of one at 0.73) | collapse chain identifiers before counting (**0% → 53%**) |
| whole COX-1 run crashed | heme outcompeted the real NSAID by atom count in ligand auto-detection | prevalence-based ligand detection — a component appearing in hundreds of PDB entries is an additive, not a ligand. Replaced an unbounded hand-curated exclusion list |
| OpenMM parser crash | `TER` records for filtered-out chains written through unconditionally | chain-scope `TER` like `ATOM`/`HETATM` |
| trypsin hold-out reported "falsified" | **my own verification error** — compared raw SIFTS numbers against classical literature numbering without checking the conventions correspond | motif-anchored calibration; result was actually confirmed |
| a literature claim contradicted the tool | an AI search summary asserted a residue assignment **not present in the paper it cited** | fetch primary sources; the tool was right |

The last two are the reason every literature check in this project quotes a primary source.

## Design notes

**Receptor preparation.** Receptor bonds are built from PDBFixer/OpenMM's `Topology` rather than
re-perceived from PDB text, because RDKit's PDB bond guessing is unreliable on real structures.
Aromaticity is assigned manually (Topology bonds carry no bond order), and catalytic metal ions are
retained explicitly (PDBFixer's `removeHeterogens` strips them).

**Why UniProt numbering internally.** Author residue numbers are not comparable across depositions.
Everything is keyed on SIFTS-mapped UniProt positions, with author numbering shown alongside in
output so results remain usable in PyMOL/ChimeraX.

## Project structure

```
posegate/
├── posegate/
│   ├── conserved_contacts.py   # mining + leave-one-out self-validation
│   ├── selectivity.py          # N-way cross-isoform comparison
│   ├── residue_mapping.py      # SIFTS -> UniProt remapping
│   ├── numbering_display.py    # multi-coordinate numbering
│   ├── autopsy.py              # per-pose interaction report + fitted score
│   ├── docking.py              # AutoDock Vina orchestration
│   └── receptor_prep.py        # Topology-based receptor building
├── scripts/
│   ├── mine_target.py                  # one-command entry point
│   ├── ca_family_matrix.py             # family-wide selectivity matrix
│   ├── baseline_comparison.py          # vs naive baselines
│   ├── selectivity_vs_experiment.py    # vs ChEMBL Ki (negative)
│   └── ligand_conditioned_test.py      # ligand-conditioned (negative)
├── PREREGISTRATION_*.md        # predictions committed before running
├── HOLDOUT_RESULT_*.md         # outcomes, including corrections
├── *_RESULT.md                 # baseline, family matrix, experimental
└── THREATS_TO_CONCLUSION.md    # what could still overturn the conclusion
```

## Testing

```bash
pytest
```

76 tests. Regression tests exist for every bug in
[Correctness history](#correctness-history) — including synthetic reproductions of the ERα
numbering failure and the CA13 chain-collapse failure.

## License

MIT. See [LICENSE](LICENSE).
