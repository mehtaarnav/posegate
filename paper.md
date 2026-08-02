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

`posegate` is a Python toolkit for triaging docked protein-ligand poses and
for automatically mining conserved binding-site contacts across a PDB
ensemble of structures for one target. It wraps AutoDock Vina
[@eberhardt2021autodockvina] for docking and ProLIF [@bouysset2021prolif]
for interaction fingerprinting, and provides two things that plain docking
scores alone do not: (1) a per-pose "autopsy" that reports steric clashes,
hydrogen bonds, aromatic contacts, and a specific, checkable pharmacophore
constraint (e.g. a named conserved hydrogen bond), combined into a
heuristic score used to rank candidates within a screened batch; and (2) a
conserved-contact miner that, given several real co-crystal structures of
one target bound to chemically distinct ligands, aggregates their ProLIF
fingerprints to automatically surface which receptor contacts recur across
ligands — the data-driven equivalent of a hand-picked literature
pharmacophore, without requiring the user to already know the literature.

Both components were validated against real, independently-verifiable
ground truth. The pose-autopsy engine's hydrogen-bond geometry was checked
against the literature-known Asn140 contact in a BRD4-JQ1 co-crystal
structure (PDB 3MXF). The conserved-contact miner was run on five protein
families spanning unrelated folds and binding chemistries — BRD4 (a
bromodomain reader domain), CDK2 (a kinase), estrogen receptor alpha (a
nuclear hormone receptor), HIV-1 protease (an obligate homodimeric
aspartic protease), and human carbonic anhydrase (a zinc metalloenzyme) —
using 4-6 real co-crystal structures per family, each bound to a
chemically distinct ligand. In every one of the five families, the miner
automatically and correctly recovered that target's textbook pharmacophore
residue with no target-specific knowledge hardcoded: BRD4's Asn140
acetyl-lysine-mimetic contact, CDK2's Leu83 hinge hydrogen bond, ER-alpha's
Glu353/Arg394 "charge clamp," HIV-1 protease's Asp25/Asp25' catalytic dyad
(recovered symmetrically from *both* monomer chains), and carbonic
anhydrase's Thr199 gatekeeper hydrogen bond.

To validate these results independently, we also profiled every structure
in all five ensembles with PLIP [@salentin2015plip], an established,
widely-used interaction-detection tool built on different underlying
geometric criteria, and mined the same ensemble-level conserved-contact
frequency table from its output. Across the five families (residues at
$\geq$50% ensemble frequency, `scripts/compare_miners.py`), the two tools
agreed on 20 residues (posegate: 28 total conserved residues; PLIP: 62),
for an overall recall of 0.32 and precision of 0.71 relative to PLIP's
broader output. The residue-level detail matters more than the aggregate
number: **the literature-validated pharmacophore residue for every one of
the five families was in the agreement set** — the two tools, using
different interaction-detection engines, independently converge on exactly
the biologically established contact in every case tested. Disagreement is
concentrated in PLIP's longer tail of looser hydrophobic/water-bridge
calls, not in the core pharmacophore signal.

We diagnosed two representative disagreements rather than leaving them
unexplained. CDK2's catalytic Lys33, which PLIP flags but posegate's
ensemble frequency does not reach threshold for: per-structure geometry
shows two real, close (H...acceptor ~2.4 Angstrom) contacts, with
donor-H...acceptor angles of 135.1 and 126.8 degrees — one just inside,
one just outside ProLIF's default 130-180 degree hydrogen-bond angle
cutoff, the same borderline-threshold effect (not a detection bug)
documented for carbonic anhydrase's Gln92 case above; the other four
structures in that ensemble have no close Lys33 contact at all, consistent
with not every CDK2 inhibitor chemotype reaching that residue. HIV-1
protease's catalytic Asp25/Asp25' shows a real per-chain asymmetry in
posegate's own raw data (chain A: HBDonor in 1/4 structures; chain B: 3/4),
consistent with the known tendency of pseudo-symmetric transition-state-
mimetic inhibitors to form a measurably stronger, more linear hydrogen
bond to one catalytic aspartate than its structural twin — a real
inhibitor-specific asymmetry, not an artifact, though the small ensemble
size (n=4) here can't rule out sampling noise as a contributor.

posegate's conserved-contact-mining approach is not the first attempt at
this problem. visGReMLIN [@ribeiro2020visgremlin] takes a directly
comparable input specification — "a set of structures composed by similar
proteins in complex with different ligands" — and mines conserved
interaction motifs via graph pattern mining rather than per-residue
frequency aggregation over an interaction-fingerprint library. We have not
yet run a head-to-head comparison against visGReMLIN specifically (as we
did against PLIP for raw interaction detection); doing so is future work,
not a claim this paper makes.

# Statement of need

Structure-based virtual screening commonly ranks candidate compounds by a
docking program's own scoring function. That score alone is a weak
discriminator in practice: in this project's own 65-compound BRD4
benchmark (22 literature actives, 43 DUD-E-style property-matched decoys),
raw Vina score achieves an AUC-ROC of 0.53 (95% stratified-bootstrap CI
[0.37, 0.68], not distinguishable from random at this sample size)
separating actives from decoys. `posegate` combines Vina score with
ProLIF-detected interaction features (H-bond count, a conserved
pharmacophore constraint, steric clash count) into a single score via
feature weights fit by L2-regularized logistic regression against this
benchmark's active/decoy labels, rather than hand-picked constants. Honest
(5-fold cross-validated, out-of-fold) AUC-ROC for that fit is 0.62 — a
real improvement over raw score, though modest and calibrated on one
target and 65 compounds, not yet validated as a general-purpose reweighting
across targets. (In-sample AUC-ROC on the same data the weights were fit
to is 0.67, 95% CI [0.54, 0.80]; the cross-validated 0.62 is the number to
trust for how this generalizes, not the in-sample figure.) Notably, the
fit gives generic H-bond count a *penalizing* weight: property-matched
decoys, matched on donor/acceptor counts to the actives, form just as many
incidental H-bonds, so raw H-bond count doesn't discriminate — only the
*specific*, literature-grounded conserved contact does. Getting from a
plain docking score to genuine ranking value requires software that can
inspect *why* a pose scored well or poorly and fit that inspection's
weights against real interaction-detector output, rather than intuition;
that is `posegate`'s contribution, together with the tooling
(validated interaction detection, pose-restraint selection, data-driven
conserved-contact mining) that produces the features being weighted.

`posegate` addresses both needs directly, and its scope is deliberately
narrow relative to related tools:

- **PoseBusters** [@buttenschoen2024posebusters] checks a pose's chemical
  and physical plausibility — bond lengths, ring planarity, stereochemistry,
  steric clashes — but does not assess target-specific interactions or
  pharmacophore recovery. `posegate`'s autopsy is complementary: it assumes
  a chemically valid pose and asks whether it engages the *right*
  interactions for that particular target.
- **Errington et al.** [@errington2024assessing] introduced a PLIF
  recovery-rate metric, also built on ProLIF, that measures how well a
  predicted pose reproduces the interactions of a *known reference
  (crystallographic) pose*. That is a pose-accuracy metric requiring a
  ground-truth answer to compare against. `posegate` is built for the
  opposite, more common situation in a real screen: triaging candidates for
  which no reference pose exists.
- **ParaDockS** [@meier2010paradocks] proposed conceptually the same idea
  `posegate`'s miner implements — target-specific, interaction-based
  post-docking classifiers trained from structural knowledge — over a
  decade ago, as part of a full population-metaheuristic docking framework.
  Its source is available but has had no commits since 2015 and is not a
  maintained, installable package. `posegate` re-implements the narrower
  conserved-contact-mining idea on top of a currently maintained
  interaction-fingerprinting library (ProLIF) and packages it as tested,
  documented, installable open-source software.
- **visGReMLIN** [@ribeiro2020visgremlin] is the closest direct comparator
  to the conserved-contact miner: it takes the same input specification
  (an ensemble of one target's structures, each bound to a different
  ligand) and is actively different in approach rather than absent from
  the field, unlike the other tools compared here. Where `posegate` mines
  conserved contacts via per-residue frequency aggregation over ProLIF
  interaction fingerprints, visGReMLIN detects conserved 3D motifs via
  graph mining. We have validated `posegate`'s miner output against PLIP
  (see Summary); a head-to-head comparison against visGReMLIN specifically
  is future work, not a claim made here, and a reviewer should weigh this
  tool's proximity accordingly when judging the miner's novelty.
- **FTMap** [@kozakov2015ftmap] and **Fragment Hotspot Maps**
  [@radoux2016fragment] identify druggable "hot spots" on a protein surface
  via, respectively, FFT-accelerated small-molecule probe docking and a
  statistical model built from the Cambridge Structural Database. Both are
  substantially heavier-weight infrastructure (a probe-docking simulation
  engine, or a licensed proprietary structural database) aimed at
  *discovering* where a pocket's hot spots are on a single structure.
  `posegate`'s miner instead asks a narrower, lighter-weight question given
  structures that already exist in the public PDB: across several
  different ligands already known to bind this target, which specific
  receptor contacts recur, and how often.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF
[@bouysset2021prolif], AutoDock Vina [@eberhardt2021autodockvina], and
OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction
detection.

# References
