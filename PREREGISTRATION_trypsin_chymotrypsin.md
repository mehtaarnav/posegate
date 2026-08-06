# Pre-registration: trypsin vs chymotrypsin selectivity comparison

Written and committed BEFORE running any structure search, mining, or
alignment for this family. The point of pre-registration is that the
prediction below cannot have been shaped by having already seen the
result -- unlike CA II/IX/XII and CDK2/CDK9, where target literature was
consulted only *after* mining, but the target *choice* and top-N list
were already in hand by the time literature was checked. This time the
prediction is committed to git history first, so it's independently
checkable that it predates the run.

## Family and reason for choosing it

Bovine trypsin (UniProt P00760) vs bovine chymotrypsin A (UniProt
P00766). Third distinct fold tested by this method: CA II/IX/XII was a
zinc hydrolase, CDK2/CDK9 was a kinase, this is a serine protease.
Chosen specifically because it has the single most famous, most
textbook-cited specificity-determinant residue in all of enzymology --
if the method can't recover this, that's a real, meaningful failure, not
a marginal one.

## Prediction (falsifiable, stated before running anything)

1. **Catalytic triad recovered as shared scaffold.** The catalytic triad
   His57/Asp102/Ser195 (chymotrypsin numbering convention) should appear
   in `shared_by_all` (or the pairwise-shared equivalent for a 2-isoform
   comparison) for both proteins -- this is the enzymatic machinery
   itself, invariant across the whole serine protease superfamily, not a
   selectivity determinant. Expected outcome, not the interesting part.

2. **The S1 pocket specificity residue should diverge and land in the
   isoform-unique or isoform-divergent bucket.** Classic biochemistry:
   trypsin has Asp189 at the base of the S1 pocket, giving it an anionic
   pocket that binds the positively charged side chains of Lys/Arg
   substrates. Chymotrypsin has a neutral/hydrophobic S1 pocket (Ser189,
   Gly216, Gly226 region) that instead prefers bulky aromatic/hydrophobic
   side chains (Phe/Tyr/Trp). This is THE textbook example of a single
   residue identity change controlling substrate specificity in an
   enzyme family. Prediction: position 189 (chymotrypsin numbering)
   should appear as a top-N conserved contact in at least one of the two
   proteins, with divergent residue identity (Asp in trypsin vs
   Ser/other in chymotrypsin) once mapped through the alignment -- the
   same shape of result as CA's Phe131/Val131 and CDK's Leu83/Cys106.

3. **Honest failure condition, stated in advance:** if position 189 (or its
   aligned equivalent) does NOT appear in either protein's top-N ligand
   contacts, or the residue identity does not diverge in the documented
   direction, this pre-registration is falsified and will be reported as
   such -- not reinterpreted after the fact to still count as a win.

## What will NOT be done

No PDB IDs have been chosen yet as of this commit. No RCSB search has
been run yet. No literature beyond general textbook knowledge already
cited above (Asp189/S1-pocket mechanism, common knowledge in any
biochemistry text, not looked up specifically for this prediction) has
been consulted. The next commit after this one will contain the actual
run and its real outcome, whatever that turns out to be.
