# Family-wide selectivity matrix: human carbonic anhydrases

Seven human CA isoforms mined independently and compared in one pass via
`scripts/ca_family_matrix.py`. This is the capability that does not
reduce to something a chemist can eyeball on a pair of structures --
see BASELINE_COMPARISON_RESULT.md for why that distinction matters (a
naive single-structure baseline recovers ~65% of pairwise selectivity
residues).

## Coverage

Of 12 catalytically relevant human CA isoforms, accessions verified
against UniProt (names and organism checked, not assumed):

| mined | excluded (too few ligand-bound structures) |
|---|---|
| CA1, CA2, CA4, CA7, CA9, CA12, CA13 | CA3 (4), CA5A (0), CA5B (0), CA6 (1), CA14 (2) |

Reliability is not uniform and is not smoothed over: CA2/CA9/CA1/CA13
are MODERATE-to-HIGH; **CA4 (n=9) and CA7 (n=8) are LOW** by the
project's own measured ensemble-size curve, and anything resting on
those two columns alone should be treated as a lead, not a finding.

## Result

8 invariant positions -- the shared catalytic scaffold, and a check on
whether the method is finding real chemistry: His94 and His119 (two of
the three zinc-coordinating histidines), His64 (the proton shuttle),
Leu197/Thr198 (gatekeeper region), Gln92, Trp5, Pro200. All eight are
textbook CA active site. None is selectivity-exploitable, and a program
that targeted them would hit every isoform at once.

10 variable positions -- candidate selectivity handles, ranked by how
many distinct amino acids occur across the family.

## Three literature-confirmed determinants recovered

| pos (this numbering) | classical | difference | status |
|---|---|---|---|
| 130 | 131 | Phe (CA2) / Val (CA9) | confirmed earlier, PMC7534198 |
| 131 | 132 | Gly (CA2) / Asp (CA9) | confirmed earlier |
| **199** | **200** | **Thr (CA2) / His (CA1)** | **confirmed here, new** |

The third is the point of this exercise. Position 199 is a known CA1/CA2
discriminator -- CA I carries a unique His200 where CA II has Thr200,
and swapping it shifts CA II's anion-inhibition Ki values toward CA I's.
**The earlier three-way CA2/CA9/CA12 comparison could not have found it,
because CA1 was not in that comparison.** Widening from three isoforms
to seven surfaced a real determinant that the narrower view structurally
could not reach.

## Position 91: checked against literature, CONFIRMED

Position 91 was flagged by this matrix as the most variable contact
position in the family -- six distinct amino acids across seven isoforms
(F/T/R/I/K/K/L) -- purely from structural mining, with no literature
input. It was recorded here as an unvalidated lead and then checked.

Verified against the primary source (Probing the Surface of Human
Carbonic Anhydrase for Clues towards the Design of Isoform Specific
Inhibitors, PMC4355338), which states:

- "Residue positions 67, 91, and 131 establish this region termed the
  *selective pocket*"
- "Residues at position 91 seem to have the highest variability, in
  terms of specific residues type and between amino acid properties
  (i.e., hydrophilicity/hydrophobicity) between isoforms"
- "Position 91 can be termed a 'hot-spot' for the design of isoform
  specific inhibitors"
- Its Table 3 lists CA II as carrying **Isoleucine** at position 91,
  matching this matrix's CA2 cell exactly.

This is a stronger form of agreement than the earlier confirmations,
which were presence/absence. Here the matrix's own ranking metric --
count of distinct amino acids across the family -- independently
reproduced the literature's characterization of position 91 as *the*
highest-variability position in the CA active site, and ranked it first
of ten without being told what to look for.

Coverage of the literature-defined selective pocket (classical
positions 67, 91, 131) is partial and stated as such: the matrix
recovered 91 and 131 (the latter as its position 130; see the numbering
note below) but **missed 67**, which was not a top-10 contact in any
isoform.

### A numbering trap worth recording

The classical-to-raw offset in this family is NOT constant. Five
independent anchors (Gln92, His94, His96, His119, Val121) confirm offset
0 in the 90-121 region, so raw 91 = classical 91. But raw 130 =
classical 131 and raw 199 = classical 200, an offset of +1 further along
-- there is an insertion between. Any comparison against CA literature
numbers has to be calibrated locally, exactly as the trypsin hold-out
required (see HOLDOUT_RESULT_trypsin_chymotrypsin.md, where assuming a
uniform offset produced a false negative).

### A source-reliability trap worth recording

The first literature search on position 91 returned a confident summary
asserting "CA I has Ile91 while CA II has Phe91" -- the reverse of what
this matrix found. Fetching the paper it cited (PMC12914371) showed that
claim is **not in that paper at all**; its only mention of position 91 is
an engineered "I91L" CA IX mimic variant, not an isoform comparison. The
matrix was right and the summary was fabricated. The direct sequence
alignment settled it independently (CA2 Ile91 aligns to CA1 Phe92, with
clean flanking matches L-L, Q-Q, F-F, H-H on both sides), and the
correct assignment was later corroborated by PMC4355338's Table 3.
Search summaries were treated as leads requiring a primary source, not
as evidence.


## Positions 134 and 19: checked against literature

**Position 134 (classical Val135): CONFIRMED, and completely.** The
literature states a design rule with explicit per-isoform residues:
"directing steric/hydrophobic bulk into the 130s subpocket (which
includes Val135 in CA II) exploits differences between isoforms, where
CA II has Phe131/Gly132/Val135 versus CA IX has Val131/Asp132/Leu135 and
CA XII has Ala131/Ser132/Ser135."

Against this matrix (raw 130/131/134 = classical 131/132/135), all nine
cells match exactly -- CA2 F/G/V, CA9 V/D/L, CA12 A/S/S. The matrix
independently reproduced a published three-residue subpocket design rule
in full, and ranked all three positions inside its top-10 variable list.

**Position 19 (classical Asp19): NOT CONFIRMED.** No source describes it
as an active-site or selectivity-relevant residue. Recorded as a false
positive. It is also the only checked candidate that was mined by a
single isoform (CA2 alone); every confirmed candidate except classical
132 was mined by three or more. That is a weak signal on n=4, not a rule,
but it is the obvious thing to test if a credibility filter is ever
added.

Running tally of candidates checked: 91 confirmed, 130 confirmed, 131
confirmed, 134 confirmed, 199 confirmed, 19 false positive. All flagged
candidates were checked and all outcomes are reported, negative included.

## Other novel candidates, still unvalidated

Positions 20, 121, 140 and 259 are variable and have not been checked
against literature.

## Limitations

- The top-10 cutoff is arbitrary; a residue at rank 11 in every isoform
  is invisible here.
- Position 130 has no aligned residue in CA4 ('-'), so that cell is a
  genuine gap, not a match.
- The catalytic zinc is a real, correctly-mined top-10 contact in CA1
  and CA13 but has no sequence position, so it cannot appear in a
  residue-identity matrix. It is reported separately. Every clinical CA
  inhibitor coordinates it, which makes it pan-isoform by definition and
  the opposite of a selectivity handle -- but its absence from the
  matrix is a limitation of the representation, not evidence it does not
  matter.
- Two of the seven columns (CA4, CA7) rest on LOW-reliability ensembles.

## A correctness bug this exercise exposed

Building this surfaced a real defect in the core miner, since fixed
(commit 90d0757): residue labels carried chain identifiers, so
independent depositions that lettered equivalent crystallographic copies
differently split one physical residue into two counted labels. CA13
scored 0% leave-one-out top-1 on a 15-structure HIGH-reliability
ensemble because PHE132.A (0.33) and PHE132.B (0.40) were counted
separately instead of as one residue at 0.73. After the fix, CA13 scores
53% and PHE132/THR200 take ranks 1-2 -- the chemically correct answer.
The pairwise comparisons never exposed this; running seven isoforms did.
