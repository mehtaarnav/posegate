# posegate — archived

**Status: development stopped 2026-08-09. Not maintained. Not recommended for use.**

A computational-chemistry toolkit that mines a target's conserved binding-site contacts from its
own co-crystal structures, self-validates by leave-one-out, and compares those contacts across the
isoforms of a protein family.

This repository is kept as a record, not as a tool. It stopped because its central claim was
tested and did not hold.

## What was established

**Supported:** the miner recovers literature-confirmed selectivity residues across several
structurally unrelated protein families, and does so better than sequence variability or naive
single-structure geometry alone.

**Tested and not supported:** that its output predicts experimentally achievable isoform
selectivity. Two pre-specified tests with controls — one against ~43,000 ChEMBL Ki measurements,
one ligand-conditioned — both failed. In the second, the control outperformed the hypothesis.

Development stopped by a decision rule fixed before those tests were run.

## Where to read

| document | contents |
|---|---|
| [`SESSION_SUMMARY.md`](SESSION_SUMMARY.md) | **start here** — full project narrative, both sessions, what was learned |
| [`EXPERIMENTAL_VALIDATION_RESULT.md`](EXPERIMENTAL_VALIDATION_RESULT.md) | the pair-level negative result and the variance decomposition explaining it |
| [`LIGAND_CONDITIONED_RESULT.md`](LIGAND_CONDITIONED_RESULT.md) | the ligand-conditioned negative result |
| [`THREATS_TO_CONCLUSION.md`](THREATS_TO_CONCLUSION.md) | audit of what could overturn the remaining claims |
| [`BASELINE_COMPARISON_RESULT.md`](BASELINE_COMPARISON_RESULT.md) | comparison against naive baselines |
| [`FAMILY_MATRIX_RESULT_carbonic_anhydrase.md`](FAMILY_MATRIX_RESULT_carbonic_anhydrase.md) | seven-isoform selectivity matrix, confirmations and one false positive |
| `PREREGISTRATION_*.md` / `HOLDOUT_RESULT_*.md` | predictions committed before running, and their outcomes |

## If you are reusing any of this

The transferable parts are methodological, not the code:

- **Verify the ensemble before trusting any number from it** — accession, organism, wild-type vs
  mutant, chain consistency, numbering convention. Most wrong answers here came from skipping one.
  An "ERα" ensemble once silently contained ERβ and ERR-γ.
- **Author residue numbers are not comparable across depositions.** Remap via SIFTS and verify the
  accession, or conserved-contact counts are meaningless.
- **Chain identifiers split one physical residue into two** when depositions letter copies
  differently. This produced a 0% accuracy score on an otherwise healthy ensemble.
- **Test against something independent of the structural record.** Literature agreement is partly
  guaranteed when the literature was written from the same structures.

`scripts/selectivity_vs_experiment.py` and `scripts/ligand_conditioned_test.py` are reusable
harnesses for checking a structural claim against experimental binding data. They are the most
portable thing here.

The code runs (76 tests, clean install verified), but it is unmaintained and its headline claim
is unsupported. Treat it as a worked example, not a dependency.

## License

MIT. See [LICENSE](LICENSE).
