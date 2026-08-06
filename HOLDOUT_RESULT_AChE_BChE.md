# Hold-out result: acetylcholinesterase vs butyrylcholinesterase

Outcome of the prediction pre-registered in `PREREGISTRATION_AChE_BChE.md`
(committed dee35e1, before any structure search or mining for this
family).

## What happened

13/14 structures prepped and mined cleanly (one, 5FPQ, hit the known
PDBFixer missing-residue-region error, already documented elsewhere in
this project -- not a new issue). One structure (6CQX) was correctly
flagged by the asymmetric-multichain warning added in commit 41015a7 --
an isolated case (1/13), not systemic like chymotrypsin's 4-6/7-12, and
did not prevent a usable result. LOO: 39% top-1, MODERATE reliability at
n=13 -- lower confidence than CA/CDK, but not degenerate the way
chymotrypsin was.

Learning directly from the trypsin numbering mistake (see
HOLDOUT_RESULT_trypsin_chymotrypsin.md's correction), the raw-vs-
classical numbering question was checked BEFORE comparing anything
against the pre-registered prediction, not after getting a apparent
mismatch. Located AChE's catalytic serine via the cholinesterase
family's conserved `GESAG` motif (classically Ser203): found at raw
position 234, giving an offset of +31 from classical numbering
(consistent with AChE's known ~31-residue signal peptide). Classical
Phe295 and Phe297 -> raw 326 and 328 -- both confirmed as Phe in the raw
sequence.

## Result against the pre-registered prediction

**Both raw positions (326, 328) are already in AChE's own top-10 mined
output** (`PHE326.A`, `PHE328.A`) -- not just present somewhere in the
full list, but empirically ranked in the top-10 by the miner without any
manual steering. Aligned to BChE via sequence alignment (Bio.Align,
BLOSUM62, same tooling as every other family comparison in this
project):

- AChE raw 326 (Phe) -> BChE raw 314 (**Leu**)
- AChE raw 328 (Phe) -> BChE raw 316 (**Val**)

This is an exact match to the pre-registered prediction on both
residues, not just one: the documented aromatic-to-aliphatic acyl-pocket
narrowing that is the textbook basis of the entire AChE/BChE selectivity
story. **Confirmed.**

## Assessment

Third independent, structurally distinct family (zinc hydrolase [CA],
kinase [CDK], serine protease [trypsin], now a different serine
hydrolase fold [cholinesterase]) where this method's core finding
replicates: the miner's own top-ranked, self-validated conserved-contact
residues land on real, independently-published selectivity-determining
positions, without being told what to look for. Unlike the trypsin case,
this one required no post-hoc correction -- the numbering-convention
check was done proactively this time, and the result matched on the
first pass.
