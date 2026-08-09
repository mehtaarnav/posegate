# posegate — project summary across both sessions

A narrative record of how this project developed, what it concluded, and
why it stopped. This is the *retrospective*; `README.md` is the tool's
public documentation and answers a different question.

Two working sessions, 57 commits, 2026-08-02 to 2026-08-09.

> **Provenance note.** Session 2 is first-hand. Session 1 ("Posegate
> project scaffold", 1836 messages) was read back from its transcript
> plus git history — direct quotes below are verbatim from it, the
> technical arc is from the commits.

---

## Session 1 — scaffold and JOSS preparation (2026-08-02 → 2026-08-05)

**How it started.** The opening request was for a *scaffold only*:

> "computational chemistry tool named 'posegate'. Do not write the
> internal algorithm logic — leave function bodies empty with simple
> `pass` or `raise NotImplementedError` so I can implement the code
> myself."

What emerged instead, over 1836 messages, was a fully implemented,
benchmarked toolkit with a JOSS paper attached. That scope drift — from
"give me the skeleton, I'll fill it in" to a submission-ready package —
is the first thing worth noticing about this project, and it set up the
pattern the second session spent its time undoing: build fast, then
discover late what the artefact could actually support.

**What was built.** A pose-triage toolkit: AutoDock Vina docking
orchestration, a per-pose "autopsy" reporting clashes, hydrogen-bond
geometry, aromatic and metal contacts, and a fitted score for ranking
candidates. Alongside it, a conserved-contact miner that derives a
target's pharmacophore from its own co-crystal structures instead of
requiring one from the literature.

**The arc.** Early commits are scaffolding and JOSS paper preparation.
Then a reframing — *"Reframe the miner as the primary component"* —
after a cross-target transferability study found that BRD4-fitted
ranking weights **do not transfer to other targets**. That was the
session's most consequential finding: the pose-ranking score, originally
the headline, was demoted to exploratory, and the miner became the
project.

**Correctness work.** A substantial fraction of the commits are fixes
found by results not looking right: decoy construction, receptor
aromaticity perception, expression-tag residues numbered below 1,
catalytic metal ions being stripped, ligands lost to embedding failures,
a cross-validation leak, stale weights. An external audit prompted a
batch of these.

**Where it ended.** Five-target validation reported, a carbonic anhydrase
inversion diagnosed, author metadata filled in, paper rebuilt against
JOSS's required structure — a project positioned for submission. Its
final exchange was the discovery that submission was impossible: the
repo went public 2026-08-02, and JOSS requires six months of prior
public development, putting earliest eligibility at ~2027-02-02. The
session closed on that, verbatim:

> "I should have caught this the first time we discussed JOSS readiness,
> not after building out the whole submission package. That's a real
> miss: I checked word counts, citation accuracy, and required paper
> sections in detail, but never checked repository-eligibility
> requirements at all."

Session 2 opened on exactly that problem.

---

## Session 2 — validation, pivot, and stop (2026-08-05 → 2026-08-09)

### JOSS is impossible, and the goal dissolves

JOSS requires 6+ months of non-concentrated public development history.
The repo had three days. Submission was not merely premature but
ineligible, and the framing that had organised the previous session
disappeared with it. What followed was a search for what the tool was
actually *for*, which took several turns and produced real casualties.

### The bug that reframed everything

Extending the miner to more structures produced **ERα leave-one-out
top-1 accuracy of 0%**. My first hypothesis — that ERα has genuinely
diverse binding modes — was wrong. Two real defects:

1. The "ERα" ensemble contained **ERβ and ERR-γ structures**, pulled in
   by keyword search and never checked against UniProt.
2. Independent depositions assign the **same author residue number to
   different physical residues**.

Fixing both — mandatory SIFTS remapping plus accession verification —
took ERα from **0% to 80%**. This set the pattern for everything after:
*the ensemble is the experiment, and an unverified ensemble produces a
confident wrong number rather than an obvious failure.*

### Pivot to selectivity

A resistance-mutation direction was tried first: mined contacts vs
Stanford HIVdb positions gave odds ratio 16.8, p = 0.0007. Then a
pocket-vs-surface control was added and it collapsed to p = 0.057. The
result was real but marginal, and the direction was dropped.

Cross-isoform **selectivity mapping** replaced it, and this is where the
project's strongest work sits.

### Pre-registration

Predictions were committed to git *before* any structure search, so the
record could not be reshaped after the fact. Two families were
pre-registered and both confirmed (AChE/BChE, COX-1/COX-2). A third
(trypsin/chymotrypsin) was first reported as **falsified** — then
corrected when the failure turned out to be **my own verification error**:
I compared raw SIFTS numbers against classical literature numbering
without checking the conventions corresponded. Motif-anchored
calibration showed the tool had been right. The wrong verdict was kept
in the file alongside the correction.

### The family matrix

Seven human carbonic anhydrase isoforms compared in one pass. All eight
invariant positions were textbook catalytic machinery. Confirmed against
primary literature: **position 91** (which the matrix ranked first by
variability, and which literature independently calls *the*
highest-variability position and a named selectivity hot-spot), the
**131/132/135 subpocket** (all nine isoform cells matching a published
design rule exactly), and **Thr200** — which the earlier three-isoform
comparison structurally could not have found, because CA1 was not in it.
One candidate, **position 19**, was checked and reported as a false
positive.

Building this exposed another correctness bug: residue labels still
carried chain identifiers, so equivalent crystallographic copies counted
as separate residues. CA13 scored **0% top-1 at HIGH reliability** —
`PHE132.A` at 0.33 and `PHE132.B` at 0.40 instead of one residue at 0.73.
Collapsing chain identifiers took it to 53%.

### Testing the claim properly

Three tests, escalating in independence:

- **vs naive baselines** — full pipeline 4/4 on the pre-registered subset
  against 65% for single-structure geometry. Underpowered (n=4, p ≈ 0.18)
  and stated as such. A useful negative: ensemble geometry (50%) is
  *worse* than single-structure geometry, because conservation ranking
  demotes the rim residues where selectivity lives.
- **vs sequence-only** — could a plain alignment find these residues with
  no structure? Position 91 ranks 11/260 by sequence variability;
  Thr200 ranks 159/260. Structure does the discriminative work. Threat
  retired.
- **vs experiment** — the decisive one.

### The negative results

The accumulated evidence could not separate *"finds
selectivity-determining positions"* from *"finds binding-site positions,
of which determinants are a subset."* Every literature confirmation was
consistent with both. So the test had to be against something
independent of the structural record.

**Pair-level**, ~43,000 ChEMBL Ki measurements, one pre-specified
statistic: ρ = +0.244, **p = 0.287**. Whole-protein control ρ = +0.199.
11.3% of random non-contact position sets matched or beat it. Post-hoc
alternative observables were *worse*, closing the escape hatch.

**Variance decomposition** explained why: **74.5% of selectivity variance
is within isoform pairs (ligand-driven)**, only 25.5% between them. The
test had collapsed each pair to a median and measured the minority
component.

**Ligand-conditioned**, the corrected version using compounds with solved
structures: main effect ρ = −0.171 (p = 0.193, wrong sign); pre-specified
control ρ = +0.341 (**p = 0.0076**). The control succeeded and the
hypothesis failed. Compounds touching *zero* divergent positions were the
most selective in the sample.

The deeper obstacle: **compounds with crystal structures are
systematically the non-selective ones.** Pan-CA sulfonamides crystallise
readily and have decades of study; genuinely selective inhibitors are
recent and largely structurally unsolved. The data needed to test the
hypothesis properly is largely absent from the PDB.

Development stopped there, by a decision rule fixed before the tests ran.

---

## Where it landed

**Supported:** posegate systematises and ranks what the structural record
already encodes — recovering literature-confirmed selectivity residues
across several protein families, better than sequence variability or
naive geometry alone.

**Tested and not supported:** that its output predicts experimentally
achievable isoform selectivity.

**Engineering state:** 76 tests, clean install verified in a fresh
environment, regression coverage for every bug below.

## The bugs, and what they have in common

| symptom | cause | fix |
|---|---|---|
| ERα 0% top-1 | wrong proteins in ensemble; incomparable author numbering | SIFTS remapping + accession check (0% → 80%) |
| CA13 0% at HIGH reliability | chain letters split one residue into two | collapse chain identifiers (0% → 53%) |
| COX-1 run crashed | heme outcompeted the real ligand by atom count | prevalence-based detection, replacing an unbounded exclusion list |
| chymotrypsin unusable | non-identical chains lettered inconsistently across depositions | detect and refuse as unsupported target class |
| HIV protease 87% → 40% | figure produced before accession verification existed | re-mined on title-verified wild-type structures |

Every one surfaced as an implausible *result*, not as a crash or a failing
test. None would have been caught by code review. The recurring lesson is
that in structure mining the ensemble is the experiment, and a bad
ensemble fails silently and confidently.

## Two errors that were mine, not the tool's

Recorded because they nearly produced wrong published conclusions:

- The **trypsin hold-out** was reported as falsified when the tool was
  right; I had compared numbering conventions without verifying they
  corresponded.
- A literature search returned a confident residue assignment that was
  **not in the paper it cited**. Fetching the primary source showed the
  tool was right and the summary was fabricated.

Both are why every literature check here quotes a primary source.

## What I would tell the next person

1. Verify the ensemble before trusting any number from it. Accession,
   organism, wild-type vs mutant, chain consistency, numbering
   convention. Most of this project's wrong answers came from skipping
   one of those.
2. Pre-register predictions. The confirmations that survive scrutiny are
   the ones committed before looking.
3. Build the test that could refute you, and run it. The most valuable
   artefact here is a reusable harness that checks structural claims
   against experimental data — and its first use returned an unfavourable
   verdict on this project's own headline claim.
4. A result that shrinks under scrutiny is the process working. The
   failure mode is shipping the version that never got scrutinised.
