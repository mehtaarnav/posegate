# Hold-out result: trypsin vs chymotrypsin

## CORRECTION (see bottom of file)

The verdict originally in this file ("does not replicate CA and CDK,"
position 189 "not confirmed") was **wrong, and the error was mine, not
the miner's**. I checked for the literal numeral "189" in the raw
SIFTS-remapped output without verifying that raw numbering actually
tracks classical chymotrypsin-numbering in this region -- it doesn't,
by a constant, motif-verifiable offset of 5. Once diagnosed properly
(see "Clean diagnosis" section below, added after the original verdict),
trypsin's own top-ranked mined residue (raw label "Asp194", 0.81
frequency, 81% LOO top-1, rank 1 of the whole ensemble) IS the real
classical-189-equivalent S1-pocket residue. The prediction holds. The
original "not confirmed" text below is kept, struck through in spirit
but not in content, so the mistake and its correction are both visible
rather than quietly edited away.

Outcome of the prediction pre-registered in
`PREREGISTRATION_trypsin_chymotrypsin.md` (committed a22f359, before any
structure search or mining for this family). Reporting the actual result,
including where it did not confirm the prediction, per that document's
own stated commitment.

## What happened

**Trypsin (P00760):** got a robust ensemble after two real, legitimate
data-quality fixes (not result-driven): 6 of the first 13 candidate
structures failed to parse in RDKit (amidine/benzamidine-group nitrogen
valence errors, a known RDKit PDB-parsing limitation on charged
amidinium groups), so the ensemble was expanded with additional verified
wild-type structures to compensate. Final result: 16/20 usable
structures, 81% LOO top-1 accuracy, HIGH reliability tier -- as robust as
the CA and CDK runs.

**Chymotrypsin (P00766):** hit a second, distinct, previously-unseen bug
during prep (`'NoneType' object has no attribute '_current_chain'`),
traced to a real defect in `prep_ensemble.py`'s chain-filtering logic
(TER records for filtered-out chains were being written through
unconditionally, crashing OpenMM's PDB parser -- fixed in commit
deaf9b2, see that commit message for detail). Fixing it raised prepping
from 4/12 to 7/8 structures. But even with that fix, mining produced an
unusable signal: 0% LOO accuracy, every top-10 residue at singleton
frequency (1/7 structures each), with the same-looking positions
(Ser195, Asp194, Gly216) scattered across chain labels A, B, and C.

## Diagnosis of the chymotrypsin failure

Chymotrypsin's mature enzyme is three disulfide-linked chains (A/B/C)
produced by proteolytic cleavage of one polypeptide -- unlike HIV
protease's homodimer (two literally identical chains, so any A/B
labeling inconsistency is harmless), chymotrypsin's chains are not
interchangeable, so a residue's chain letter is part of its identity.
The mined data is consistent with different PDB depositions assigning
chain letters to chymotrypsin's three fragments inconsistently (e.g. the
catalytic Ser195 appearing as both `SER195.A` and `SER195.C` across
different structures), which fragments the same physical residue's
identity across the ensemble exactly the way the ERalpha numbering bug
fragmented residue identity across depositions -- except this is a
chain-level version of that problem, not a residue-number version, and
the current pipeline has no mechanism to detect or correct for it. This
is a real, distinct limitation, not fixed as part of this task.

## Verdict against the pre-registered prediction

- **Catalytic triad as shared scaffold**: untestable. Chymotrypsin's
  ensemble never produced a usable signal to compare against trypsin's.
- **S1-pocket specificity residue (Asp189, chymotrypsin numbering)
  recovered as a divergent top-N contact**: NOT confirmed on the side
  that could be checked. Trypsin's full mined output (not just top-10)
  was checked across the entire 185-200 region; position 189 does not
  appear at any frequency in any structure. The pre-registration's own
  stated failure condition applies here directly.

**Original (WRONG) conclusion, kept for transparency:** "This hold-out
test does not replicate the CA and CDK results... reported as a genuine
partial failure, not reinterpreted post-hoc as a win."

## Clean diagnosis (correction)

Checked whether trypsin's raw SIFTS-remapped numbering actually tracks
the classical chymotrypsin-numbering convention used in the pre-
registration's prediction, rather than assuming it does. It does not,
directly: trypsin's UniProt canonical sequence (P00760) begins
`MKTFIFLALLGAAVAFPVDDDDKIVGGYTC...` -- a signal peptide plus the
trypsinogen activation propeptide (the classic `DDDDK` enteropeptidase
cleavage site) before the mature enzyme's own N-terminus (`IVGG...`,
the famous "Ile16" of trypsin activation biochemistry) even begins. Raw
UniProt position 189 is nowhere near the real classical position 189.

A constant-offset correction doesn't work either (serine protease
classical numbering has non-linear insertion-code loops), so the
correct method is to anchor on an unambiguous, motif-identified
landmark instead of arithmetic: the catalytic serine's `GDSGGP` motif,
universal across the whole serine protease superfamily. Located
precisely (`seq.find('GDSGGP')`, catalytic Ser at the motif's 3rd
position) in both proteins' raw sequences:

- Trypsin: catalytic Ser at raw 200. Six residues N-terminal (classical
  195 - 189 = 6): raw 194 = **Asp**.
- Chymotrypsin: catalytic Ser at raw 195 (coincidentally close to its
  own classical number here, unlike trypsin). Six residues N-terminal:
  raw 189 = **Ser**.

Both match the pre-registered prediction exactly: anionic Asp in
trypsin's S1 pocket, neutral Ser in chymotrypsin's. And empirically,
raw-194 IS trypsin's own miner output's rank-1 residue by a wide margin
(0.81 frequency, 13/16 structures, 81% LOO top-1) -- the exact numeric
label ("194") that made this look like a miss when compared naively
against "189" was itself the artifact, not the underlying result.

**Corrected verdict:** the S1-pocket-divergence prediction IS confirmed
on the trypsin side, both by independent sequence anchoring and by the
miner's own robust, high-confidence output. The catalytic-triad-shared-
scaffold half of the prediction remains untestable (chymotrypsin's own
ensemble is still unusable -- that diagnosis was correct and stands;
see below). This hold-out is a genuine additional replication of the
method's core finding, not a failure -- once evaluated correctly rather
than checked with a numeric-matching shortcut. The lesson that survives
from the original (wrong) analysis: don't compare raw SIFTS-remapped
labels against literature numbers without verifying the numbering
conventions actually correspond -- exactly the class of error this whole
project exists to catch when it happens in the *miner's* output, and it
turned out to happen in my own verification code instead.

## What remains genuinely unresolved

Chymotrypsin's ensemble is still unusable (0% LOO, singleton
frequencies) due to the diagnosed chain-labeling inconsistency across
its three non-symmetric chains -- that finding was not affected by the
numbering correction above and stands as a real, current limitation.
See "Diagnosis of the chymotrypsin failure" above for detail.
