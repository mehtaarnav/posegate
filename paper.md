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

`posegate` is a Python toolkit whose primary component is a conserved-contact miner. Given an
ensemble of one target's co-crystal structures, each bound to a chemically distinct ligand, it
aggregates ProLIF interaction fingerprints [@bouysset2021prolif] and reports which receptor
contacts recur across the ensemble and at what frequency. This yields a target-specific
pharmacophore constraint for targets where the literature does not already supply one. The miner
was evaluated on five protein families with unrelated folds, cross-checked against PLIP
[@salentin2015plip], and compared against the published results of visGReMLIN
[@ribeiro2020visgremlin], the closest existing tool.

A second component applies the mined constraint to docked poses. It wraps AutoDock Vina
[@eberhardt2021autodockvina], reports the clashes, hydrogen bonds and aromatic contacts each pose
makes and whether it satisfies the conserved contact, and combines these into a score fitted
against labeled data. That score's ranking performance is modest and does not transfer between
targets; we report it as a characterization of the approach's limits rather than as a validated
screening method.

# Statement of need

Interaction-aware analysis of docked poses requires a target-specific pharmacophore to check
against: a named receptor contact that genuine binders are expected to make. For well-studied
targets this comes from the literature, as we initially took BRD4's Asn140 contact. For a
less-studied target it is often unavailable, and identifying one from the structural literature
is precisely the expertise a newcomer to that target lacks. Many such targets nonetheless have
several ligand-bound structures already in the PDB, which collectively contain the answer.
`posegate`'s miner extracts it, requiring no prior knowledge of the target beyond a list of its
PDB entries and their bound ligands.

The downstream motivation is that raw docking score is a weak discriminator. On this project's
90-compound BRD4 benchmark (30 literature actives, 60 DUD-E-style property-matched decoys), raw
Vina score separates actives from decoys with an AUC-ROC of 0.60 (95% stratified-bootstrap CI
[0.49, 0.72]), not clearly distinguishable from random at this sample size. Whether a pose makes a
specific, mechanistically meaningful contact is a more direct question, and the pose-triage
component checks exactly that. The intended audience is researchers working on a target with no
established literature pharmacophore, and those running structure-based screens who want
per-pose interaction detail behind a docking score.

# State of the field

Rescoring docked poses by the interactions they make, rather than by the docking score alone, is
an established approach. Structural interaction fingerprints [@deng2004sift] encode a pose's
contacts as a bit string for comparison against known complexes, and SPLIF [@da2014splif] uses
that comparison to recover actives that the docking score alone would reject. `posegate`'s
restraint-guided pose selection is the same idea applied at selection time rather than after it:
among the poses Vina returns, prefer one that satisfies a required contact. What differs is where
the required contact comes from. These methods take a reference complex, or a pharmacophore
supplied by the user, as given; `posegate` derives it by mining the target's own deposited
structures, which is what makes it applicable to a target whose reference pharmacophore has not
been established.

`posegate` occupies a narrow scope relative to several other tools. PoseBusters
[@buttenschoen2024posebusters] checks a pose's chemical and physical plausibility but not
target-specific interaction recovery, so the two are complementary. Errington et al.
[@errington2024assessing] introduced a ProLIF-based metric for how closely a predicted pose
reproduces a known reference pose, which requires ground truth, whereas `posegate` addresses the
screening case where none exists. ParaDockS [@meier2010paradocks] proposed the same idea the
miner implements over a decade ago, but its source has had no commits since 2015 and is not an
installable package. FTMap [@kozakov2015ftmap] and Fragment Hotspot Maps [@radoux2016fragment]
locate hot spots on a single structure, by probe docking and by a model built from the Cambridge
Structural Database, requiring much heavier infrastructure for a different question.

visGReMLIN [@ribeiro2020visgremlin] is the closest direct comparator, taking the same input but
mining conserved 3D motifs by graph pattern mining rather than by per-residue frequency
aggregation. It was released only as a web server, without source code, and both advertised URLs
are currently unreachable, so we compared against its published CDK case study instead
(`scripts/compare_visgremlin.py`). Of the 26 binding-site atoms of Schonbrunn et al., spanning 9
residues, visGReMLIN recovered 18 atoms across 8 residues; on a 22-structure CDK2 ensemble
`posegate`'s miner reports contacts at all 9 residues, or 8 of 9 counting only
non-van-der-Waals interactions, missing the same residue visGReMLIN missed, HIS84. The hinge
residue PHE82, which visGReMLIN recovers through an aromatic motif, is recovered here as a
hydrophobic contact rather than as pi-stacking. The scores are not directly comparable, since
visGReMLIN's is atom-level over 73 complexes. The comparison also showed that ensemble
composition matters more than size: adding 14 structures from one fragment-screen series removed
ASP145 from the output entirely. The README reports both results in full.

# Software design

Interaction detection is centralized on ProLIF rather than implemented directly. Receptor bond
perception proved to be a correctness risk rather than a formality: RDKit's native PDB parser
supplements template-based bonding with a distance-based heuristic that can introduce spurious
bonds at tight turns, and with that heuristic disabled its residue-template matcher can leave
most of a chain unbonded without reporting an error, as we observed on an unremarkable structure
where only the first two residues of 298 were bonded. `posegate.receptor_prep` therefore builds
the receptor molecule from PDBFixer/OpenMM's `Topology.bonds()` rather than inferring
connectivity from PDB text. Docking adds a ligand-size-aware search box and restraint-guided pose
selection, since Vina has no restraint term; the README documents both.

# Research impact statement

`posegate` has not yet been used in a published research study. In all five protein families
tested, the conserved-contact miner recovered the established literature pharmacophore with no
target-specific knowledge hardcoded: BRD4's Asn140 acetyl-lysine-mimetic contact, CDK2's Leu83
hinge hydrogen bond, estrogen receptor alpha's Glu353/Arg394 charge clamp, HIV-1 protease's
Asp25/Asp25' catalytic dyad, recovered symmetrically from both chains of this obligate homodimer,
and carbonic anhydrase's Thr199 gatekeeper hydrogen bond. Cross-checking against PLIP
[@salentin2015plip] placed the literature-validated residue in the agreement set in every family.
Agreement across all conserved residues was recall 0.32 and precision 0.71 against PLIP's broader
output, with disagreement concentrated in PLIP's tail of looser calls; two representative
disagreements trace to a ProLIF angle cutoff and to a real per-structure asymmetry rather than to
detection errors, and the README gives both diagnoses in full.

We fitted the pose-triage score independently on each of the five families, each with its own
DUD-E-style property-matched benchmark (64 to 261 compounds) and its own mined pharmacophore
constraint substituted for the feature that constraint checks. Raw-Vina baselines ranged from
0.52 (CDK2) to 0.72 (estrogen receptor alpha), consistent with the target-dependent variance
reported for docking generally [@trannguyen2020litpcba] rather than indicating a defect in any
one setup. Fitted weights are not comparable in magnitude across targets at these sample sizes,
so we compare the direction of each feature's association with activity instead, fitting each
target's weights on standardized features (the scaler refit on each cross-validation fold's own
training data rather than on the full dataset beforehand, so no test-fold statistics leak into
training) and bootstrapped over 200 resamples per target, reported only when the sign is stable
in at least 90% of them (`scripts/compare_feature_weights.py`). The mined conserved contact is
the one feature whose
direction holds across all five targets, with bootstrap sign stability between 0.94 and 1.00.
Generic hydrogen-bond count and the raw docking score itself both reverse direction between
targets. This is consistent with the miner, rather than the fitted score, being the component
of `posegate` that generalizes across targets, and is the strongest evidence for that available
from this study.

Carbonic anhydrase is a useful illustration of why. Its raw-Vina baseline was 0.25, meaning Vina
ranked actives *worse* than decoys on average, yet its mined Thr199 contact still discriminated
correctly (bootstrap-stable at 1.00). We traced the inversion rather than treating it as noise:
docked poses place carbonic anhydrase inhibitors 1.4-2.9 Angstrom from the catalytic zinc,
consistent with genuine metal coordination, against 1.8-4.2 Angstrom for decoys, so pose
selection is finding the correct site. The zinc ion carries zero partial charge in the prepared
receptor, and AutoDock Vina's own documentation states that it disregards atomic charges on
metal ions during scoring; the AutoDock4Zn extension exists specifically to address this
[@santosmartins2014autodock4zn], and we did not use it here. Vina's raw score for this target is
therefore unreliable by a documented limitation of the scoring function, not a defect in our
receptor preparation, and the conserved-contact feature discriminates correctly regardless,
because it does not depend on that score.

Two further limitations bound how far the pose-triage results should be read. Property-matched
decoy sets in the style of DUD-E [@mysinger2012dude] carry analogue and decoy bias that a fitted
model can learn in place of learning protein-ligand interaction: Chen et al.
[@chen2019hiddenbias] traced deep-learning enrichment on DUD-E to exactly that artifact rather
than to generalization. Any margin our fitted scores show over their raw-docking baselines
should therefore be treated as an upper bound, and confirming it would require a benchmark built
to avoid these biases, such as LIT-PCBA [@trannguyen2020litpcba]. This caveat applies to the
fitted weights and not to the miner, which uses no decoys at all.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF [@bouysset2021prolif], AutoDock Vina
[@eberhardt2021autodockvina], and OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction detection.

# References
