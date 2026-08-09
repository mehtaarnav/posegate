# What could still overturn the conclusion

The conclusion under test: *PoseGate is a useful conserved-contact and
selectivity-mapping tool.* Convenience, docs, packaging, scalability and
automation are deliberately out of scope here -- only things that could
change whether the claim is true.

## Retired by direct test

**Could sequence variability alone find these residues, with no
structure?** This would have been fatal: the matrix ranks candidates by
count of distinct amino acids across the family, which is a *sequence*
property. If a plain alignment reproduced the answer, the entire
structural pipeline would be decoration.

Tested. Per-position variability across the same seven isoforms, whole
protein, no structural input:

| position | classical | sequence-only rank |
|---|---|---|
| 91 | 91 | **11 / 260** |
| 131 | 132 | 46 / 260 |
| 130 | 131 | 90 / 260 |
| 134 | 135 | 91 / 260 |
| 199 | 200 | **159 / 260** |

Sequence alone buries the confirmed determinants. Thr200 -- a verified
CA1/CA2 discriminator -- ranks 159th. The most sequence-variable
positions protein-wide (24, 85, 182, 9, 34, 36, 37, 55, 58, 60...) are
overwhelmingly surface loops; exactly one known determinant (91) appears
in the sequence-only top 15. The structural matrix's top 10 contains at
least four confirmed determinants.

The structural filter is doing the discriminative work: it reduces 260
positions to ~18 genuine contact positions, and variability ranking is
only meaningful *within* that set. **Threat retired.**

## Tested and NOT survived

**Does the tool's output predict experimentally achievable selectivity?**
This is the direct H1-vs-H2 test the whole threat list turned on, and it
was built and run (EXPERIMENTAL_VALIDATION_RESULT.md, ~43,000 ChEMBL Ki
measurements, 21 isoform pairs). Mined contact divergence vs experimental
selectivity: rho = +0.244, p = 0.287. Whole-protein divergence control:
rho = +0.199. Random non-contact position sets match or beat the real one
11.3% of the time. Alternative observables are worse, not better, so the
null is not an artifact of the pre-specified choice.

At the aggregate level the evidence does **not** favour H1 over H2.
Threat D is therefore no longer a framing concern to be argued about --
it is a measured result. The supportable claim is "systematises and ranks
what the structural record encodes", not "predicts achievable
selectivity".

## (1) Conclusion-threatening

**A. Retrospective validation with an unknown denominator.** Most
confirmations followed run -> inspect output -> pick a residue -> check
literature. Even position 91, which was recorded as a lead before being
checked, was picked *because* it ranked first. The danger is a garden of
forking paths: if failures go unrecorded, a string of confirmations
means little. Partially mitigated now -- of the three candidates flagged
in the matrix writeup, all three were checked and one (position 19,
classical Asp19) is reported as a false positive. But the discipline is
young and fragile, and the historical record before it was adopted has
an unknown denominator.

**B. The margin over naive baselines is not established.** The
pre-registered comparison is 4/4 for the full pipeline against ~65% for
single-structure geometry, n=4, p ~ 0.18. Every claim of the form
"better than doing the obvious dumb thing" rests on this and it is not
significant. This is the most directly fixable threat: more
pre-registered families.

**C. Never tested on the actual use case.** Every validation family was
selected *because* it has literature to check against. The claimed value
is for targets that lack exactly that. This is close to
unfalsifiable-by-construction: you cannot validate on targets defined by
having no ground truth. Partially addressable by simulating the
condition -- pre-register predictions on a family whose literature is
deliberately not consulted until after -- but never fully.

**D. Ensembles encode medicinal-chemistry intent, not just pocket
biology.** Deposited structures exist because someone designed ligands
for that site. "Conserved across the ensemble" therefore partly measures
what a design program aimed at. Decades of sulfonamide chemistry against
CA means the mined contacts partly reflect sulfonamide design decisions.
This does not threaten the claim of *systematizing what is known*; it
does threaten any claim of *discovery independent of prior knowledge*,
and the two have been used somewhat interchangeably.

## (2) Worth noting, non-threatening

- **Reliability tiers derived from one target.** The ensemble-size curve
  came from CDK2 alone and is applied to every family.
- **LOO measures internal consistency, not correctness.** A homogeneous
  ensemble (one med-chem series) scores high and means little. Observed
  in practice on a CDK2 fragment series.
- **Top-N cutoff is arbitrary.** Position 67, a genuine member of the
  literature's three-residue selective pocket, was missed at top-10.
- **False-positive rate is unquantified.** One known FP (position 19) out
  of four candidates checked. Four is not a rate.
- **Prevalence threshold (200) was set from ~8 observed compounds.**
- **Scope is enzymes only.** No GPCRs, ion channels, or protein-protein
  interfaces.
- **Two of seven CA columns rest on LOW-reliability ensembles** (CA4 n=9,
  CA7 n=8).
- **The tool can emit contradictory signals** -- CA13 reported HIGH
  reliability and 0% top-1 simultaneously, which required human
  interpretation to resolve as a bug rather than a finding.

## (3) Irrelevant to the conclusion

- Asymmetric multi-chain receptors unsupported -- detected, refused, and
  disclosed. A bounded scope limit, not a correctness threat.
- Multi-residue and peptide ligand detection -- same.
- Numbering-convention friction -- an interpretation cost, now surfaced
  directly in the output.
- obabel / PDBFixer install friction -- operational.
- The stale Chen et al. benchmark -- superseded by stronger evidence.

## Honest summary

Threat A is methodological and can be fixed by discipline already
started. Threat B is statistical and fixable with more pre-registered
families. Threat D is a framing problem: the evidence supports
"systematizes and ranks what the structural record already encodes"
much more strongly than it supports "discovers what nobody knew."
Threat C may be permanent.


## Post-mortem: the ensemble-QC pivot, checked and closed

After development stopped, one pivot looked strong enough to check: drop
selectivity, keep the machinery that WAS validated, and reframe as a
structural-ensemble QC tool -- verify that a set of PDB structures are
the same protein, consistently numbered, correctly ligand-assigned and
internally coherent, before anyone trains a model or runs ensemble
docking on them. Every one of those failures bit this project, on
extremely well-studied proteins.

Checked before committing to it. It is occupied:

- **PDBCleanV2** (bioRxiv 2025.02.14.638326) "compares the sequences of
  all chains in each structure to a set of reference sequences" and
  "standardizes chain names and numbering", explicitly rectifying
  "mislabeling errors from original structure submission to the PDB".
  That is the ERalpha identity contamination and the CA13 chain-letter
  split, both handled -- by MUSCLE alignment against reference sequences
  rather than SIFTS lookup, but to the same effect.
- **LP-PDBBind** and **PDBbind CleanSplit** (Nat Mach Intell 2025)
  already address structure-dataset contamination and train/test leakage
  for the ML use case.
- **wwPDB OneDep** covers per-structure validation.

The only remaining differentiator would be the ground-truth-free
leave-one-out coherence signal, which is narrow and only meaningful for
contact mining specifically. Not a product.

That makes five ideas checked across this project -- JOSS eligibility,
the curated-database coverage gap, resistance prediction, selectivity
prediction, ensemble QC -- and five that came back occupied, confounded
or falsified. Four of the five could have been checked in under an hour
BEFORE building. That is the single most expensive lesson here.
