# Hold-out result: COX-1 vs COX-2

Outcome of the prediction pre-registered in `PREREGISTRATION_COX1_COX2.md`
(committed dee35e1, before any structure search or mining for this
family).

## Deviations from the pre-registration, disclosed as they happened

1. **No human COX-1 structures exist in the PDB at all** (0 hits).
   Nearly all classic COX-1 structures, including the ones that
   established the Ile523/Val523 mechanism historically, are ovine
   (P05979). Switched species before any mining -- a data-availability
   constraint discovered at the very first step, not a result-driven
   choice.
2. **Chain symmetry check** (promised as the first execution step):
   the built-in asymmetric-multichain detector (commit 41015a7) flagged
   1/7 COX-2 structures (5F1A) in isolation, not systemically -- same
   pattern as AChE's isolated flag, not chymotrypsin's systemic one.
   Did not block a usable result.
3. **Two new ligand-detection bugs found and fixed along the way**,
   same discipline as every prior hold-out: COX is a heme-dependent
   peroxidase (HEM/COH/MNH outcompeted real inhibitors by atom count)
   and a membrane-associated, glycosylated protein (BOG detergent, NAG
   glycosylation sugar did the same). Both fixed in LIKELY_NON_LIGAND
   with tests before re-running, not worked around by hand-picking
   structures.

## Result against the pre-registered prediction

COX-1 (ovine, 12/13 usable, 83% LOO top-1, MODERATE reliability):
`ILE523.A` appears directly in the top-10 mined output (0.33 frequency,
4/12 structures) -- matching the pre-registered raw number exactly, no
offset needed this time.

COX-2 (human, 7/7 usable, 71% LOO top-1, LOW reliability given small n):
raw position 523 is Asn, not Val -- but COX-2 is 4 residues longer than
COX-1 (604 vs 600), so raw numbers aren't directly comparable across the
two proteins, the same lesson as every prior cross-protein comparison in
this project. Aligned COX-1's confirmed raw-523 (Ile) through
Bio.Align/BLOSUM62 to COX-2's own numbering: **raw 509 (Val)** --
already present in COX-2's own top-10 mined output (`VAL509.A`, rank 10,
0.14 frequency, 1/7 structures -- a weak signal given the small ensemble,
but present, not absent).

**Confirmed on both sides**, matching the pre-registered prediction
exactly: Ile in COX-1, Val in COX-2, at the correctly aligned position --
the textbook mechanistic basis for the entire selective COX-2 inhibitor
drug class (celecoxib, rofecoxib, etc.), independently recovered by the
miner without being told what to look for.

## Assessment

Fifth structurally distinct family (zinc hydrolase, kinase, serine
protease, cholinesterase, now a heme-dependent membrane-associated
peroxidase) where the core finding replicates. The COX-2 side is weaker
evidence than the others (n=7, LOW reliability, the target residue at
rank 10 not rank 1) -- worth stating plainly rather than folding into
the same confidence as COX-1's cleaner result. A larger COX-2 ensemble
would be the natural next check if this family's result needed to bear
more weight than "one more confirmation among several."
