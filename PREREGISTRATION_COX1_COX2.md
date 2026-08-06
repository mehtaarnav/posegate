# Pre-registration: COX-1 vs COX-2

Written and committed BEFORE running any structure search, mining, or
alignment for this family.

## Family and reason for choosing it

Human cyclooxygenase-1 (COX-1/PTGS1, UniProt P23219) vs cyclooxygenase-2
(COX-2/PTGS2, UniProt P35354). Chosen as the second of two hold-out
families for this round, alongside AChE/BChE (see
PREREGISTRATION_AChE_BChE.md), specifically because it is arguably the
single most famous structure-based drug-selectivity story in modern
pharmacology (the rational design basis for celecoxib, rofecoxib, and
the whole "coxib" NSAID class), giving an unusually high-confidence,
independently-checkable prediction. COX exists as a homodimer; before
treating this as inside the supported domain, the dimer's chain
symmetry will be checked empirically (via the same
is_asymmetric_multichain check added in commit 41015a7) as the very
first step of execution, before any interpretation of mined results --
if it turns out non-symmetric, this family is out of the supported
domain and will be reported as such rather than forced through.

## Prediction (falsifiable, stated before running anything)

COX-1 and COX-2 have near-identical catalytic active sites, but COX-2 has
a single amino acid substitution -- Ile523 in COX-1 versus the smaller
Val523 in COX-2 -- that opens a "side pocket" absent in COX-1. This
single-residue difference is the textbook basis for the entire selective
COX-2 inhibitor drug class.

Prediction: position 523 (COX numbering) should appear as a top-N
conserved ligand contact in at least one of the two proteins' own mined
output, with the aligned position showing divergent residue identity
(Ile in COX-1, Val in COX-2) once mapped through sequence alignment.

## Honest failure condition, stated in advance

If position 523 (or its properly motif/alignment-verified equivalent --
not assumed via raw number, per the lesson from the trypsin correction)
does not appear in either protein's top-N mined contacts, or the
residue identity does not diverge in the documented direction, this is
falsified and will be reported as such.

## What will NOT be done before the next commit

No PDB IDs chosen yet. No RCSB search run yet. No structure-specific
literature beyond the general textbook mechanism cited above consulted.
