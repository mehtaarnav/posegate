# Baseline comparison: does the pipeline beat something dumber?

Every prior validation in this project showed the miner RECOVERS known
selectivity residues. None showed the machinery was necessary -- the
confirmed residues are binding-site residues, and "list what's near the
ligand" returns binding-site residues too. This tests that directly.

Run via `scripts/baseline_comparison.py` against the nine already-prepped
ensembles from this project's five validated families.

## Methods compared

| | interaction typing | ensemble |
|---|---|---|
| **B1** geometry, single structure | no | no |
| **B2** ProLIF, single structure | yes | no |
| **B3** geometry, ensemble | no | yes |
| **FULL** the actual pipeline | yes | yes |

Metric: recovery@10 -- does the literature-confirmed selectivity residue
appear in the method's top 10? B1/B2 are single-structure methods, so
each is run on every ensemble member and reported as a mean rate rather
than depending on an arbitrary choice of which structure.

## Results

```
Target                B1 geom/1  B2 prolif/1  B3 geom/ens    FULL
CA II (P00918)             58%          33%           no     YES
CA IX (Q16790)             73%          18%          YES     YES
CDK2 (P24941)             100%          87%          YES     YES
CDK9 (P50750)              92%          85%          YES     YES
Trypsin (P00760)           50%          94%          YES     YES  [prereg]
AChE (P22303)              23%          23%          YES     YES  [prereg]
COX-1 (P05979)            100%          33%           no     YES  [prereg]
COX-2 (P35354)             86%          14%           no     YES  [prereg]
MEAN (all 8)               73%          48%          62%    100%
MEAN (prereg only)         65%          41%          50%    100%
```

## The selection-bias correction

The "all 8" row overstates FULL. Four targets (CA II, CA IX, CDK2, CDK9)
were identified by reading FULL's own top-10 output and only then
verified against literature -- FULL cannot miss them, by construction.
Only the four marked `[prereg]` were predicted from published literature
BEFORE mining ran on that family (see `PREREGISTRATION_*.md`). The
`MEAN (prereg only)` row is the one that measures anything.

The conclusion survives the correction: on the four honest targets, FULL
recovers 4/4 while the best baseline recovers ~65%.

## What this does and does not establish

**Does:** the full pipeline outperformed every baseline on both subsets,
and the two components are complementary rather than redundant --
interaction typing alone (B2, 41%) and ensemble conservation alone (B3,
50%) each underperform their combination (FULL, 100%). That is a
coherent mechanistic result, not just a score.

**A specific, useful negative finding:** ensemble geometry (B3, 50%) is
*worse* than single-structure geometry (B1, 65%). Ranking by "fraction
of structures where this residue is near the ligand" favors the deep,
always-occupied pocket core and pushes rim residues down -- and
selectivity-determining residues are frequently at the rim. Naive
conservation actively harms this task. FULL avoids the trap because
ProLIF weights by interaction quality rather than mere proximity, so a
rim residue forming real interactions across chemically diverse ligands
still ranks high.

**Does NOT:** establish this with statistical confidence. The honest
subset is n=4. Observing 4/4 when the best baseline's rate is 65% has
p ~ 0.18 under a binomial null -- directionally consistent, nowhere near
significant. Four more pre-registered families would be needed before
this margin could be called established rather than suggestive.

**Also worth stating plainly:** B1, the dumbest possible method -- one
structure, raw distance, no interaction typing, no ensemble -- gets the
answer about two-thirds of the time. The pipeline's advantage is real
but it is an advantage over a fairly effective naive approach, not over
nothing. Anyone choosing whether to use this should weigh it against
that, not against zero.
