# Hold-out result: trypsin vs chymotrypsin

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

**This hold-out test does not replicate the CA and CDK results.** Per
the pre-registration's commitment, this is reported as a genuine
partial failure, not reinterpreted post-hoc as a win. What it did
produce, honestly: two real infrastructure bugs found and fixed (the
RDKit amidine-parsing sensitivity, worth knowing about for any future
serine-protease-family work, and the TER-record chain-filtering crash,
now fixed and tested), and clear evidence that the method's two prior
clean successes do not trivially generalize to every protein family --
specifically, multi-chain, non-symmetric proteins with inconsistent
cross-deposition chain labeling are a real, current blind spot, not
covered by the SIFTS residue-numbering fix that solved the analogous
problem for single-chain proteins (ERalpha) or symmetric homodimers (HIV
protease).
