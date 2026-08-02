---
title: 'posegate: ProLIF-based pose triage and PDB-ensemble conserved-contact mining for structure-based virtual screening'
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

`posegate` is a Python toolkit for two related tasks in structure-based drug discovery: judging
whether a docked protein-ligand pose is trustworthy, and discovering which receptor contacts a
target's real binders consistently engage. Docking programs like AutoDock Vina
[@eberhardt2021autodockvina] output a pose and a numeric score, but the score alone says little
about *why* a pose scored the way it did. `posegate` wraps Vina with
[ProLIF](https://github.com/chemosim-lab/ProLIF)-based interaction fingerprinting
[@bouysset2021prolif] to produce a per-pose "autopsy" — steric clashes, hydrogen bonds, aromatic
contacts, and a specific, checkable pharmacophore constraint — combined into a score fit against
labeled data rather than hand-picked. Separately, its conserved-contact miner takes an ensemble
of a target's real co-crystal structures, each bound to a chemically distinct ligand, and
automatically surfaces which receptor contacts recur across them: the data-driven equivalent of
a hand-picked literature pharmacophore, without requiring the user to already know the
literature. Both components were validated against real, independently-verifiable ground truth
across five protein families with unrelated folds, and cross-checked against an established,
independently-implemented interaction detector (PLIP [@salentin2015plip]).

# Statement of need

Raw docking score is a weak discriminator in practice. In this project's own 65-compound BRD4
benchmark (22 literature actives, 43 DUD-E-style property-matched decoys), raw Vina score
achieves an AUC-ROC of 0.53 (95% stratified-bootstrap CI [0.37, 0.68]) separating actives from
decoys — not distinguishable from random at this sample size. Getting from a plain docking score
to genuine triage value requires software that can inspect *why* a pose scored well or poorly,
using real interaction-detector output rather than the score alone. A second, related problem is
that this kind of interaction-aware triage typically needs a target-specific pharmacophore to
check against, and that pharmacophore is usually either hand-picked from the literature (as we
first did for BRD4's Asn140 contact) or simply unavailable for a less-studied target. `posegate`
addresses both: its autopsy module produces interaction-aware, fitted pose scores instead of raw
docking scores, and its conserved-contact miner derives the pharmacophore constraint itself from
existing PDB structures rather than requiring the user to already know it. The target audience is
researchers running structure-based virtual screens who want more than a docking score to rank
candidates, and researchers working on a target without an established literature pharmacophore
who want one mined directly from that target's own structural history.

# State of the field

`posegate`'s scope is deliberately narrow relative to several related tools. PoseBusters
[@buttenschoen2024posebusters] checks a pose's chemical and physical plausibility — bond lengths,
ring planarity, stereochemistry, clashes — but not target-specific interaction recovery;
`posegate`'s autopsy is complementary, assuming a chemically valid pose and asking whether it
engages the *right* interactions for that target. Errington et al.
[@errington2024assessing] introduced a PLIF recovery-rate metric, also built on ProLIF, that
measures how well a predicted pose reproduces a *known reference crystallographic pose*'s
interactions — a pose-accuracy metric requiring ground truth to compare against. `posegate` is
built for the opposite, more common screening situation: triaging candidates for which no
reference pose exists. ParaDockS [@meier2010paradocks] proposed conceptually the same idea
`posegate`'s miner implements — target-specific, interaction-based post-docking classifiers
trained from structural knowledge — over a decade ago, as part of a full population-metaheuristic
docking framework; its source is available but has had no commits since 2015 and is not a
maintained, installable package. visGReMLIN [@ribeiro2020visgremlin] is the closest direct
comparator to the conserved-contact miner, taking the same input specification (an ensemble of
one target's structures, each bound to a different ligand) but mining conserved 3D motifs via
graph pattern mining rather than per-residue frequency aggregation over interaction fingerprints;
we have validated `posegate`'s miner output against PLIP (see Summary and Research impact,
below) but have not yet run a head-to-head comparison against visGReMLIN specifically, and a
reviewer should weigh this tool's proximity accordingly when judging the miner's novelty. FTMap
[@kozakov2015ftmap] and Fragment Hotspot Maps [@radoux2016fragment] identify druggable "hot
spots" on a protein surface via, respectively, FFT-accelerated small-molecule probe docking and a
statistical model built from the Cambridge Structural Database — both substantially
heavier-weight infrastructure aimed at *discovering* where a pocket's hot spots are on a single
structure, versus `posegate`'s narrower, lighter-weight question given structures that already
exist in the public PDB: across several different ligands already known to bind this target,
which specific receptor contacts recur, and how often.

# Software design

Interaction detection is centralized on ProLIF rather than hand-rolled geometry: an earlier
implementation computed hydrogen-bond and aromatic-contact geometry directly, and migrating to
ProLIF both simplified that code and gave access to a wider, independently-maintained set of
interaction definitions. A less obvious design choice concerns receptor bond perception, which
turned out to be a real correctness risk rather than a formality: RDKit's native PDB parser
either supplements template-based bonding with a distance-based heuristic that can invent
spurious bonds on tight turns, or, with that heuristic disabled, can silently leave the large
majority of a multi-residue chain completely unbonded (observed directly: only the first one or
two residues of a 298-residue chain bonded, the rest orphaned, on an unremarkable structure).
`posegate.receptor_prep` avoids re-guessing connectivity from written PDB text altogether by
building the receptor's RDKit molecule directly from PDBFixer/OpenMM's own `Topology.bonds()` —
bonds OpenMM already computes correctly as part of running molecular simulations. Docking
combines a ligand-size-aware search box (sized to each ligand's own conformer extent plus a fixed
margin, since a fixed box much smaller than a given ligand causes catastrophic clashes) with
restraint-guided pose selection: Vina has no restraint term in its scoring function, so rather
than accepting only its single top-ranked pose, `posegate` requests several candidate poses and
selects the best-scoring one that actually satisfies a required contact, falling back to the
top-ranked pose if none do. Finally, because a fixed absolute score threshold only makes sense
relative to whatever distribution a given receptor and scoring setup actually produces, batch
screening ranks candidates by percentile within the screened set rather than a fixed cutoff.

# Research impact statement

`posegate` has not yet been used in a published research study; it is presented here as a
validated, ready-to-use tool rather than as one already integrated into an ongoing project. The
substance of its validation is what an interested user or reviewer should weigh: the
conserved-contact miner correctly and automatically recovered the textbook literature
pharmacophore, with no target-specific knowledge hardcoded, in all five protein families tested
— BRD4's Asn140 acetyl-lysine-mimetic contact, CDK2's Leu83 hinge hydrogen bond, estrogen
receptor alpha's Glu353/Arg394 "charge clamp," HIV-1 protease's Asp25/Asp25' catalytic dyad
(recovered symmetrically from both monomer chains of this obligate homodimer), and carbonic
anhydrase's Thr199 gatekeeper hydrogen bond. Cross-checking against PLIP [@salentin2015plip], an
independently-implemented interaction detector, across all five ensembles found that the
literature-validated pharmacophore residue was in the agreement set in every one of the five
families, and that overall agreement across all conserved residues was recall 0.32 / precision
0.71 relative to PLIP's broader output, with disagreement concentrated in PLIP's longer tail of
looser calls rather than in the core pharmacophore signal. Two representative disagreements were
diagnosed rather than left unexplained: a missed CDK2 lysine contact traced to a real hydrogen-bond
donor-H...acceptor angle of 126.8 degrees, just outside ProLIF's default 130-180 degree cutoff on
an otherwise genuine close contact; and an HIV-1 protease catalytic-aspartate chain asymmetry
traced to a real per-structure asymmetry in the underlying data, consistent with known
pseudo-symmetric-inhibitor binding behavior rather than a detection artifact. The screening/pose-
ranking component of the tool is weaker and reported as such rather than polished: its fitted
score reaches a cross-validated AUC-ROC of 0.62 against a raw-Vina baseline of 0.53 on one
65-compound benchmark, a real but modest improvement not yet shown to generalize across targets.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF [@bouysset2021prolif], AutoDock Vina
[@eberhardt2021autodockvina], and OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction detection.

# References
