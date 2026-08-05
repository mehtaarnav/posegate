---
title: 'posegate: mining conserved protein-ligand contacts from PDB ensembles, with interaction-aware pose triage'
tags:
  - Python
  - cheminformatics
  - molecular docking
  - drug discovery
  - structural biology
authors:
  - name: Arnav Mehta
    orcid: 0009-0008-1688-6906
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 03 August 2026
bibliography: paper.bib
---

# Summary

`posegate` is a Python toolkit for structure-based drug discovery. Its main component mines a
target's own deposited protein structures to find which receptor contacts its known binders
share, without requiring any prior knowledge of that target. Its second component checks docked
poses against that constraint and ranks them accordingly.

Structure-based screening asks whether a candidate molecule's docked pose looks like a real
binder's. Answering that well requires a pharmacophore: a specific receptor contact, such as a
hydrogen bond to a named residue, that genuine binders are known to make. For a well-studied
target this pharmacophore comes from published literature. For a less-studied target it usually
does not exist yet, and working it out by hand requires exactly the domain expertise a researcher
new to that target lacks. `posegate`'s miner solves this by extracting the pharmacophore directly
from data: given several structures of one target, each bound to a different ligand, it reports
which receptor contacts recur across them, using [ProLIF](https://github.com/chemosim-lab/ProLIF)
[@bouysset2021prolif] to detect interactions and PDBFixer/OpenMM to prepare each receptor. A
second component applies the resulting constraint to poses from
[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) [@eberhardt2021autodockvina],
reporting whether a pose satisfies it and combining that with other pose features into a fitted
ranking score.

# Statement of need

Researchers running a structure-based screen on a target with no established literature
pharmacophore currently have two options: read the structural literature closely enough to
identify one by hand, or screen without one and rely on docking score alone. The first requires
expertise specific to that target; the second is weak. On a 90-compound BRD4 benchmark built for
this project, raw Vina score separated known actives from property-matched decoys with an AUC-ROC
of 0.60 (95% bootstrap CI [0.49, 0.72]) — not clearly better than chance at this sample size.
`posegate` gives such a researcher a third option: point the miner at the target's own PDB
entries and derive the pharmacophore from the structures that already exist. The intended
audience is researchers screening a target without a known pharmacophore, and researchers who
want per-pose interaction detail behind a docking score rather than the score alone.

# State of the field

PoseBusters [@buttenschoen2024posebusters] checks whether a pose is chemically and physically
plausible — bond lengths, ring geometry, clashes — but not whether it engages the right target
interactions; the two tools answer different questions and can be used together. Errington et al.
[@errington2024assessing] measure how closely a predicted pose reproduces a known reference pose's
interactions, which requires that reference pose as ground truth. `posegate` addresses the more
common case where no reference pose exists. ParaDockS [@meier2010paradocks] proposed the same
underlying idea `posegate`'s miner implements, over a decade ago, but its source has had no
commits since 2015 and is not an installable package.

visGReMLIN [@ribeiro2020visgremlin] is the closest existing tool: it takes the same input, an
ensemble of one target's structures each bound to a different ligand, and identifies conserved
motifs by graph mining rather than by the frequency-based approach used here. It was released
only as a web server, and neither of its advertised URLs currently resolves, so we could not run
it directly; we instead compared against the results it published for its own CDK case study. Of
26 reference binding-site atoms spanning 9 residues, visGReMLIN reported 18 atoms across 8
residues; `posegate`'s miner, applied to a 22-structure CDK2 ensemble, reports all 9 residues, and
8 of 9 when restricted to specific (non-van-der-Waals) contacts — missing the same residue
visGReMLIN missed. The two scores are not directly comparable, since visGReMLIN's is atom-level
over 73 complexes and `posegate` has no atom-level output, but the result indicates the
lighter-weight method reaches the same conclusion. FTMap [@kozakov2015ftmap] and Fragment Hotspot
Maps [@radoux2016fragment] locate druggable hot spots on a single structure via probe docking or
a statistical model, respectively; both require substantially heavier infrastructure to answer a
different question than the one `posegate` asks of structures that already exist in the PDB.

# Software design

Interaction detection is centralized on ProLIF rather than computed by hand, which was the
project's original approach; this both simplified the code and gave access to a wider,
independently maintained set of interaction definitions. Receptor preparation turned out to carry
real correctness risk. RDKit's native PDB bond perception either invents spurious bonds at tight
turns or, with that heuristic disabled, silently fails to bond most of a multi-residue chain.
`posegate.receptor_prep` avoids both failure modes by building the receptor molecule directly from
PDBFixer/OpenMM's own computed bonds. Those bonds carry no order, so nothing in a prepared
receptor is aromatic unless assigned explicitly; `receptor_prep` assigns it from residue templates
for the standard aromatic side chains, which is what allows pi-stacking to be detected at all.
Catalytic metal ions are retained rather than stripped, since a metalloenzyme's ligand may bind by
coordinating one directly; because no interaction-fingerprinting library models a coordination
bond as such, that contact is instead reported geometrically, by distance, alongside the
fingerprint-derived ones.

Docking wraps AutoDock Vina with a search box sized to each ligand's own extent, since a fixed box
across a chemically diverse library causes clashes for larger ligands. Vina's scoring function has
no restraint term, so `posegate` requests several candidate poses per ligand and selects the
best-scoring one that satisfies the mined constraint, falling back to the top-ranked pose when
none does. Ligands are docked concurrently across a process pool, since Vina's own internal
threading scales sublinearly and independent worker processes make better use of many cores.

# Research impact statement

`posegate` has not yet been used in a published study; the evidence here is from validating it
directly. The miner was tested on five protein families with unrelated folds and different
pharmacophore chemistries — BRD4, CDK2, estrogen receptor alpha, HIV-1 protease, and carbonic
anhydrase — and recovered the literature-established pharmacophore in all five, with no
target-specific knowledge built into the code. An independent interaction detector, PLIP
[@salentin2015plip], placed that same literature residue in agreement with `posegate` in every
family.

The pose-ranking component is weaker evidence and is reported as such. Fitted on each target's own
benchmark, its cross-validated AUC-ROC improves on the raw-Vina baseline for four of the five
targets and is flat on the fifth (CDK2). Raw-Vina baselines themselves range from 0.25 to 0.72
across targets, a spread consistent with docking performance being known to vary sharply by
target rather than indicating a defect in any one setup. Because confidence intervals at these
sample sizes overlap too much for a pairwise AUC comparison across targets to mean anything, we
instead compare the *direction* of each feature's association with activity, bootstrapped over
200 resamples per target: the mined conserved contact is the only feature whose direction holds
across all five targets, while generic hydrogen-bond count and the raw Vina score itself both
reverse between targets. Carbonic anhydrase makes the case concretely. Its raw-Vina baseline is
inverted — actives score worse than decoys on average — because Vina's scoring function has no
term for metal coordination and its zinc-binding inhibitors depend on exactly that; docked poses
still place those inhibitors correctly, 1.4–2.9 Å from the catalytic zinc, and the mined contact
still discriminates correctly despite the unreliable raw score. Taken together, this is evidence
that the miner's output, not the fitted score, is the part of `posegate` that generalizes across
targets.

Two caveats bound these results. The benchmarks use property-matched decoys in the style of
DUD-E [@mysinger2012dude], which are known to carry bias a fitted model can learn in place of
genuine interaction signal [@chen2019hiddenbias]; the fitted-score margins above should be read as
an upper bound pending validation on a bias-corrected benchmark such as LIT-PCBA
[@trannguyen2020litpcba]. This caveat does not apply to the miner, which uses no decoys.

# AI usage disclosure

Generative AI (Claude, Anthropic) was used throughout this project's development, under direct
human direction and with human review and correction at each step. The author made all target,
scope, and methodological decisions: which protein families to validate against, which
comparators to include, what counts as sufficient evidence for a claim, and where a result should
be reported as exploratory rather than established. The AI executed implementation and debugging
under that direction — including the receptor-preparation fixes described above, the benchmark
and docking pipeline, and diagnosis of several defects found during validation (among them, an
aromaticity-detection gap that made pi-stacking unreportable on any structure, and a
differential-dropout bug in ligand preparation that was silently biasing benchmark composition) —
and drafted most of the software's documentation and this paper's prose, revised repeatedly under
author direction. All reported results were produced by running the project's own code and are
reproducible from the repository; the final paper and codebase were reviewed and verified by the
author before submission.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF [@bouysset2021prolif], AutoDock Vina
[@eberhardt2021autodockvina], and OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction detection.

# References
