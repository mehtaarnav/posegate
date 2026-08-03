---
title: 'posegate: mining conserved protein-ligand contacts from PDB ensembles, with interaction-aware pose triage'
tags:
  - Python
  - cheminformatics
  - molecular docking
  - drug discovery
  - structural biology
authors:
  - name: TODO Author Name
    orcid: TODO
    affiliation: 1
affiliations:
  - name: TODO Affiliation
    index: 1
date: TODO
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
65-compound BRD4 benchmark (22 literature actives, 43 DUD-E-style property-matched decoys), raw
Vina score separates actives from decoys with an AUC-ROC of 0.53 (95% stratified-bootstrap CI
[0.37, 0.68]), not distinguishable from random at this sample size. Whether a pose makes a
specific, mechanistically meaningful contact is a more direct question, and the pose-triage
component checks exactly that. The intended audience is researchers working on a target with no
established literature pharmacophore, and those running structure-based screens who want
per-pose interaction detail behind a docking score.

# State of the field

`posegate` occupies a narrow scope relative to several related tools. PoseBusters
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

The pose-triage component is evaluated on two targets rather than five. Its fitted score reaches
a cross-validated AUC-ROC of 0.62 against a raw-Vina baseline of 0.53 on the BRD4 benchmark it
was fitted on. Applied unmodified to an equivalent 51-compound CDK2 benchmark, with CDK2's mined
Leu83 contact substituted for Asn140, it gives 0.35 against a raw-Vina baseline of 0.37, so the
weights do not transfer. The cause is a sign reversal in one feature: generic hydrogen-bond count
is penalized on BRD4, where property-matched decoys form more incidental hydrogen bonds than
actives, but is a genuine positive signal on CDK2, where actives form more. The mined conserved
contact and the clash count keep their direction on both targets. Two targets cannot establish
which feature types are target-general and we claim no more than that, but the mined constraint
being among the features that hold is consistent with the miner, rather than the fitted score,
being the component that generalizes. Equivalent benchmarks on the remaining three families are
the natural next step.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF [@bouysset2021prolif], AutoDock Vina
[@eberhardt2021autodockvina], and OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction detection.

# References
