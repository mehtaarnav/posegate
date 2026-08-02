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
whether a docked protein-ligand pose is trustworthy, and identifying which receptor contacts a
target's known binders consistently engage. Docking programs such as AutoDock Vina
[@eberhardt2021autodockvina] return a pose and a numeric score, but that score gives little
indication of which structural features produced it. `posegate` wraps Vina with interaction
fingerprinting from [ProLIF](https://github.com/chemosim-lab/ProLIF) [@bouysset2021prolif] to
produce a per-pose report covering steric clashes, hydrogen bonds, aromatic contacts and a named
pharmacophore constraint, and combines these into a score whose weights are fitted against
labeled data rather than chosen by hand. Its conserved-contact miner takes an ensemble of a
target's co-crystal structures, each bound to a chemically distinct ligand, and reports which
receptor contacts recur across them, which provides a pharmacophore constraint for targets where
one has not already been established in the literature. Both components were evaluated against
independently verifiable ground truth across five protein families with unrelated folds, and
cross-checked against PLIP [@salentin2015plip], a separately implemented interaction detector.

# Statement of need

Raw docking score is a weak discriminator in practice. On this project's 65-compound BRD4
benchmark (22 literature actives, 43 DUD-E-style property-matched decoys), raw Vina score
separates actives from decoys with an AUC-ROC of 0.53 (95% stratified-bootstrap CI [0.37, 0.68]),
which is not distinguishable from random at this sample size. Improving on it requires software
that reports the interactions underlying a pose rather than the score alone. Interaction-aware
triage of this kind in turn needs a target-specific pharmacophore to check against, and such a
pharmacophore is typically either taken from the literature, as we initially did for BRD4's
Asn140 contact, or unavailable for a less studied target. `posegate` addresses both requirements:
its autopsy module produces fitted, interaction-aware pose scores, and its conserved-contact
miner derives the pharmacophore constraint from existing PDB structures. It is intended for
researchers running structure-based virtual screens who need more than a docking score to rank
candidates, and for those working on targets with no established literature pharmacophore.

# State of the field

`posegate` occupies a narrow scope relative to several related tools. PoseBusters
[@buttenschoen2024posebusters] checks a pose's chemical and physical plausibility, covering bond
lengths, ring planarity, stereochemistry and clashes, but not target-specific interaction
recovery. `posegate`'s autopsy module assumes a chemically valid pose and asks instead which
interactions it forms with the target, so the two are complementary. Errington et al.
[@errington2024assessing] introduced a PLIF recovery-rate metric, also built on ProLIF, that
measures how closely a predicted pose reproduces the interactions of a known reference
crystallographic pose. That metric requires ground truth for comparison, whereas `posegate` is
built for the more common screening case in which no reference pose exists. ParaDockS
[@meier2010paradocks] proposed the same idea `posegate`'s miner implements, namely
target-specific interaction-based post-docking classifiers trained from structural knowledge,
over a decade ago and as part of a full population-metaheuristic docking framework. Its source is
available but has received no commits since 2015 and is not a maintained, installable package.
visGReMLIN [@ribeiro2020visgremlin] is the closest direct comparator to the conserved-contact
miner, taking the same input specification (an ensemble of one target's structures, each bound to
a different ligand) but mining conserved 3D motifs by graph pattern mining rather than by
per-residue frequency aggregation over interaction fingerprints.
visGReMLIN was released only as a web server, without source code or a distributable package, and
both advertised URLs are currently unreachable (`vagner.dti.ufv.br/visgremlin4` refuses
connections and `homepages.dcc.ufmg.br/~alexandrefassio/gremlin/` returns PHP errors), so we could
not run it on our own ensembles. We instead compared against its published results on its own CDK
case study (`scripts/compare_visgremlin.py`). That case study scores motif recovery against the
experimentally determined CDK binding site of Schonbrunn et al., comprising 26 atoms across 9
residues; visGReMLIN recovered 18 of the 26 atoms (69%), distributed over 8 of the 9 residues.
Applied to a 22-structure CDK2 ensemble, `posegate`'s miner reports contacts at all 9 reference
residues, including HIS84, which visGReMLIN's motifs did not recover. Counting only specific
(non-van-der-Waals) interactions, it reports 7 of 9, omitting HIS84 and the hinge residue PHE82.
PHE82 is a genuine limitation rather than a scoring artifact: visGReMLIN identified its aromatic
contacts, whereas ProLIF registers the same contacts only as van der Waals proximity. The two
scores are not directly comparable, since visGReMLIN's is an atom-level score over 73 complexes
and `posegate` produces no atom-level output. We therefore claim only that the residue-level
method recovers the same published binding site at lower computational cost, not that it resolves
finer structure.

Repeating the comparison over ensembles of 6, 19 and 22 structures exposes a limitation of
frequency-based mining. Enlarging the ensemble from 6 to 19 structures degraded the result, with
ASP145 disappearing from the output altogether. Fourteen of the added structures come from a
single fragment-screen deposition series that occupies only the hinge subpocket, and restricting
that series to its drug-like members did not restore ASP145, which indicates that chemotype
homogeneity rather than ligand size is responsible. Residue coverage was comparatively stable
across all three ensembles (8 to 9 of 9), but the frequencies assigned to peripheral residues were
not. Those frequencies describe the chemistry of the ensemble rather than the pocket itself, so
ensembles for this method should be curated for scaffold diversity rather than simply enlarged.

FTMap [@kozakov2015ftmap] and Fragment Hotspot Maps [@radoux2016fragment] identify druggable hot
spots on a protein surface, using FFT-accelerated small-molecule probe docking and a statistical
model built from the Cambridge Structural Database respectively. Both require substantially
heavier infrastructure and address a different question, namely where the hot spots of a pocket
lie on a single structure. `posegate`'s miner instead operates on structures already deposited in
the public PDB and reports which receptor contacts recur, and at what frequency, across several
ligands already known to bind the target.

# Software design

Interaction detection is centralized on ProLIF rather than implemented directly. An earlier
version of the code computed hydrogen-bond and aromatic-contact geometry itself; moving to ProLIF
simplified that code and provided a wider, independently maintained set of interaction
definitions. Receptor bond perception proved to be a correctness risk rather than a formality.
RDKit's native PDB parser supplements template-based bonding with a distance-based heuristic that
can introduce spurious bonds at tight turns, and with that heuristic disabled its residue-template
matcher can leave most of a multi-residue chain unbonded without reporting an error. We observed
the latter on an unremarkable structure, where only the first one or two residues of a
298-residue chain were bonded and the remainder were left as isolated atoms.
`posegate.receptor_prep` therefore does not infer connectivity from PDB text at all, and builds
the receptor's RDKit molecule from PDBFixer/OpenMM's `Topology.bonds()`, which OpenMM computes in
order to run molecular simulations.

Docking uses a ligand-size-aware search box, sized to each ligand's conformer extent plus a fixed
margin, because a box substantially smaller than the ligand produces severe clashes. Pose
selection is restraint-guided: Vina's scoring function has no restraint term, so `posegate`
requests several candidate poses and selects the best-scoring pose that satisfies a required
contact, falling back to the top-ranked pose when none does. For batch screening, candidates are
ranked by percentile within the screened set rather than against a fixed absolute cutoff, since
the score distribution depends on the receptor and scoring setup.

# Research impact statement

`posegate` has not yet been used in a published research study, and is presented here as a
validated tool rather than one already embedded in an ongoing project. In all five protein
families tested, the conserved-contact miner recovered the established literature pharmacophore
without target-specific knowledge hardcoded: BRD4's Asn140 acetyl-lysine-mimetic contact, CDK2's
Leu83 hinge hydrogen bond, estrogen receptor alpha's Glu353/Arg394 charge clamp, HIV-1 protease's
Asp25/Asp25' catalytic dyad, recovered symmetrically from both monomer chains of this obligate
homodimer, and carbonic anhydrase's Thr199 gatekeeper hydrogen bond. Cross-checking against PLIP
[@salentin2015plip] across all five ensembles placed the literature-validated pharmacophore
residue in the agreement set in every family. Agreement across all conserved residues was recall
0.32 and precision 0.71 relative to PLIP's broader output, with disagreement concentrated in
PLIP's longer tail of looser calls rather than in the pharmacophore residues themselves. We traced
two representative disagreements to their causes. A missed CDK2 lysine contact corresponds to a
hydrogen-bond donor-H...acceptor angle of 126.8 degrees, just outside ProLIF's default 130-180
degree cutoff on an otherwise genuine close contact. An apparent chain asymmetry in the HIV-1
protease catalytic aspartates reflects a per-structure asymmetry in the underlying data, which is
consistent with known pseudo-symmetric-inhibitor binding behavior rather than a detection
artifact. The pose-ranking component is weaker: its fitted score reaches a cross-validated AUC-ROC
of 0.62 against a raw-Vina baseline of 0.53 on the 65-compound BRD4 benchmark it was fitted on, a
modest improvement.

We then tested whether the BRD4-fitted weights transfer to a second target, constructing an
equivalent 51-compound CDK2 benchmark (17 ChEMBL actives, 34 property-matched decoys) and
applying the same weights unmodified, with CDK2's Leu83 hinge contact, itself identified by the
miner, substituted for BRD4's Asn140. The weights do not transfer. AUC-ROC is 0.35 (95% CI [0.20,
0.51]), below the raw-Vina baseline of 0.37 on this benchmark. Raw Vina score is itself
sub-random on CDK2 here, so part of the shortfall is attributable to this docking run rather than
to the weights, but the weights measurably worsen it, for an identifiable reason. On BRD4, decoys
property-matched on donor and acceptor count formed more incidental hydrogen bonds than actives,
so the fit assigns hydrogen-bond count a penalizing weight. On CDK2 the same benchmark
construction produces the opposite pattern, with actives averaging 2.0 hydrogen bonds against 1.4
for decoys, so the BRD4-derived penalty acts against a signal that is genuine on this target. The
conserved-contact-hit and clash-count features retain the same relationship to activity on both
targets, whereas generic hydrogen-bond count does not. Whether a given structural feature
predicts activity is thus itself target-dependent, and weights fitted on one target's docking
output should not be assumed to transfer. Feature weights need to be refitted per target, or
restricted to features such as the conserved-contact hit whose relationship to activity is
mechanistically general rather than empirically fitted.

# Acknowledgements

This project depends on RDKit [@landrum2016rdkit], ProLIF [@bouysset2021prolif], AutoDock Vina
[@eberhardt2021autodockvina], and OpenMM/PDBFixer [@eastman2017openmm]. Validation used PLIP
[@salentin2015plip] as an independent cross-check on interaction detection.

# References
