# Ligand-conditioned test — also negative

Run via `scripts/ligand_conditioned_test.py`. This was the decisive
experiment: the pair-level test failed, the variance decomposition showed
74.5% of selectivity variance is ligand-driven, and this conditions on
the ligand instead of averaging over it.

## Design (both statistics pre-specified before running)

For each compound with a crystal structure in the mined ensemble, we know
exactly which positions that compound touches. For every isoform pair
where it also has measured Ki:

- **Hypothesis** — positions contacted that DIFFER between the two
  isoforms should predict |ΔpKi|. A compound touching nothing divergent
  has no structural basis to discriminate them.
- **Control** — positions contacted that are IDENTICAL should NOT
  predict selectivity. If it does, the signal is a contact-count size
  artifact, not a selectivity mechanism.

## Result

| | ρ | p |
|---|---|---|
| **Main** — divergent contacts | **−0.171** | 0.193 |
| **Control** — invariant contacts | **+0.341** | **0.0076** |

n = 60 observations, 11 distinct compounds.

**The hypothesis failed and the control succeeded** — the worst available
outcome. The main effect is not significant *and points the wrong way*.

Median |ΔpKi| by divergent-contact count:

| divergent contacts | n | median \|ΔpKi\| |
|---|---|---|
| 0 | 41 | **1.29** |
| 1 | 15 | 0.87 |
| 2 | 1 | 0.58 |
| 3 | 3 | 1.13 |

Compounds touching **no** divergent position are the *most* selective in
this dataset. That directly contradicts the mechanism.

The significant control is best read as a size artifact: most positions
are conserved, so invariant-contact count is close to total contact
count, and "larger, more buried compounds are more selective" is
ordinary medicinal chemistry unrelated to the tool's output.

## Limitations, stated but not used as rescue

- **Only 11 distinct compounds.** n=60 observations are heavily
  non-independent; effective n is nearer 11. Both p-values, including
  the significant control, are unreliable.
- **27/66 ligand codes mapped to ChEMBL**, and the mapped ones skew
  toward classic drugs (acetazolamide, methazolamide, ethoxzolamide) —
  compounds specifically designed as pan-CA inhibitors. Testing a
  selectivity mechanism on compounds engineered not to be selective is a
  weak test.
- **The predictor is nearly constant**: 41/60 observations have zero
  divergent contacts.

These weaken the result but do not reverse it. A pre-specified test was
run and did not support the hypothesis, while its own control did.

## The deeper problem this exposes

The compounds that have crystal structures are systematically the
non-selective ones. Classic pan-CA sulfonamides crystallise readily and
have been studied for decades; genuinely isoform-selective inhibitors are
recent, few, and largely unsolved structurally. So the data required to
test this hypothesis properly is largely absent from the PDB, and no
amount of methodological care fixes that from the current inputs.

## Verdict

Both the pair-level and the ligand-conditioned tests are negative. The
selectivity-prediction thesis is not supported, and the pre-committed
decision rule was that a negative here means stop rather than iterate.

What survives unchanged: the tool recovers literature-confirmed
selectivity residues across families, does so better than sequence
variability or naive geometry alone, and is correct engineering. What is
not supported, at either level of analysis: that its output predicts
experimentally achievable selectivity.
