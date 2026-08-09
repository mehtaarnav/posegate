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

## Novel candidates, unvalidated

Position **91** is the most variable contact position in the family --
six distinct amino acids across seven isoforms (F/T/R/I/K/K/L), mined as
a top-10 contact by CA1, CA13 and CA4. Nothing in this project has
checked it against literature. It is the obvious first candidate if this
analysis were to be pushed further, and it is stated here as a lead, not
a result.

Positions 19, 20, 134, 121, 140 and 259 are likewise variable and
unvalidated.

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
