# Pre-registration: acetylcholinesterase vs butyrylcholinesterase

Written and committed BEFORE running any structure search, mining, or
alignment for this family.

## Family and reason for choosing it

Human acetylcholinesterase (AChE, UniProt P22303) vs human
butyrylcholinesterase (BChE, UniProt P06276). Both are monomeric,
single-chain, well-crystallized cholinesterases -- chosen specifically
to stay inside the pipeline's confirmed-supported domain (single chain
near the ligand) after diagnosing chymotrypsin's asymmetric-multichain
failure mode (see HOLDOUT_RESULT_trypsin_chymotrypsin.md and commit
41015a7). Verifying this domain-fit before writing the prediction is a
target-selection step, not a peek at results -- no structures have been
fetched and no mining has run as of this commit.

## Prediction (falsifiable, stated before running anything)

The classic AChE/BChE selectivity mechanism, textbook-level established:
AChE's acyl-binding pocket is narrowed by two aromatic residues, Phe295
and Phe297, which BChE lacks -- BChE has smaller aliphatic residues
(Leu286 and Val288) at the aligned positions, giving it a larger, more
permissive acyl pocket. This is why BChE can hydrolyze bulkier
substrates (like butyrylcholine) that AChE cannot, and is a major target
of selective inhibitor design.

Prediction: at least one of these two positions (AChE numbering 295/297,
BChE numbering 286/288) should appear as a top-N conserved ligand
contact in AChE's own mined output, with the aligned BChE position
showing divergent (smaller/aliphatic) residue identity once mapped
through sequence alignment -- the same shape of result as CA's
Phe130/Val, CDK's Leu83/Cys106, and trypsin's Asp194/Ser (see
conversation and HOLDOUT_RESULT_trypsin_chymotrypsin.md for the trypsin
case, including the correction of my own earlier numbering-comparison
error there).

## Honest failure condition, stated in advance

If neither AChE position (295 or 297) appears in AChE's own top-N mined
contacts, or the aligned BChE positions do not show smaller/aliphatic
residue identity, this is falsified and will be reported as such -- not
reinterpreted after the fact, and not "rescued" by a numbering-
convention excuse without first verifying via the same motif-anchoring
discipline used to correct the trypsin result (i.e. any claimed
numbering mismatch must be independently checked via sequence, not
assumed).

## What will NOT be done before the next commit

No PDB IDs chosen yet. No RCSB search run yet. No structure-specific
literature beyond the general textbook mechanism cited above (common
pharmacology knowledge, not looked up specifically for this prediction)
consulted.
