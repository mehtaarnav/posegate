# First test against independent experimental data — negative

Run via `scripts/selectivity_vs_experiment.py` against ~43,000 ChEMBL Ki
measurements across seven human CA isoforms.

## The bottleneck this was built to remove

Every validation before this one asked *"is known residue X in the
tool's top-N?"* — binary, one residue at a time, cherry-pickable, and
answered against literature written by people who looked at the same
crystal structures the tool mines. Agreement was partly guaranteed by
construction.

So the accumulated evidence could not separate:

- **H1** — the tool identifies selectivity-*determining* positions
- **H2** — the tool identifies binding-site positions, and selectivity
  determinants are a subset of binding-site positions

Under H2 the tool is much less valuable, and every confirmation obtained
so far — position 91, the 131/132/135 subpocket, Thr200 — is equally
consistent with H2. Nothing in the previous design could tell them apart.

## The test

Experimental Ki is an independent standard: measured in wet labs, tens of
thousands of compound/isoform pairs, not dependent on anyone's structural
interpretation. If the mined positions determine selectivity, then isoform
pairs differing at more of them should show larger experimental selectivity
spreads. One pre-specified statistic over all 21 pairs at once.

## Result

| | Spearman ρ | p |
|---|---|---|
| **Main** — mined contact divergence | **+0.244** | **0.287** |
| C1 — whole-protein divergence | +0.199 | 0.387 |
| C2 — 2000 random non-contact sets | null mean +0.057 | **11.3% match or beat** |

**The prediction failed.** The direction is right and the mined positions
edge out both controls, but nothing approaches significance, and roughly
one in nine random non-contact position sets performs as well.

### The null is robust, not an artifact of the observable

Post-hoc check of alternative observables (exploratory, explicitly not
confirmatory):

| observable | ρ | p |
|---|---|---|
| median \|ΔpKi\| — *pre-specified* | +0.244 | 0.287 |
| 90th percentile \|ΔpKi\| | +0.143 | 0.536 |
| max \|ΔpKi\| | **−0.100** | 0.667 |

The pre-specified observable was the *most favourable* of the three. The
result cannot be rescued by choosing a different one.

### A real power limitation, stated but not used as an excuse

Contact divergence takes only 6 distinct values (range 3–8) across 21
pairs — a heavily compressed, heavily tied predictor with limited power.
This is a genuine weakness of the operationalisation. It does not
convert a null into a positive.

## What this does and does not mean

**Does not mean** the confirmed residues are wrong. Positions 91,
131/132/135 and 200 rest on independent structural and mutational work,
and that evidence is unaffected.

**Does mean** that *counting how many mined contact positions differ
between two isoforms does not predict how much selectivity is
experimentally achievable between them.* At this level of aggregation,
the evidence does not favour H1 over H2.

Plausible reasons, offered as interpretation and not as defence:
selectivity is unlikely to be additive in the number of differing
positions (one Phe→Val can outweigh five conservative swaps); achievable
selectivity also tracks medicinal-chemistry effort per pair, which is
enormous for CA2/CA9 and near-zero for CA13; and most measured compounds
are non-selective zinc-binding sulfonamides.

## Standing conclusion, revised

The defensible claim is narrower than before:

> PoseGate systematises and ranks what the structural record already
> encodes — recovering literature-confirmed selectivity residues, and
> doing so better than sequence variability or naive geometry alone.

The claim it does **not** support:

> PoseGate's output predicts experimentally achievable isoform
> selectivity.

That second claim was never explicitly made, but it is the one a reader
would naturally infer from "selectivity-mapping tool," and it is now
tested and unsupported.

## What the machinery is worth

The test is reusable and pre-specifiable. Any future selectivity claim
can now be checked against experiment *before* being believed, rather
than against literature derived from the same structures. That capability
did not exist before this, and its first use returned an unfavourable
answer about the project's own headline claim — which is the point of
building it.
