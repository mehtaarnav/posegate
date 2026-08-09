# posegate/scripts/selectivity_vs_experiment.py
"""Tests the tool's selectivity claim against experimental binding data
instead of against literature.

THE BOTTLENECK THIS REMOVES
---------------------------
Every validation this project has done asks "is known residue X in the
tool's top-N?" That question is binary, one residue at a time,
cherry-pickable, and -- worst -- it is answered against literature
written by people who looked at the same crystal structures the tool
mines. Agreement is therefore partly guaranteed by construction.

Consequence: the accumulated evidence cannot distinguish

  H1  the tool identifies selectivity-DETERMINING positions
  H2  the tool identifies binding-site positions, and selectivity
      determinants happen to be a subset of binding-site positions

Under H2 the tool is far less valuable, and every confirmation obtained
so far is equally consistent with H2. Nothing in the current design can
tell them apart, because there is no way to score the WHOLE output
against a standard independent of the structural record.

THE TEST
--------
Experimental inhibition constants are such a standard. They are measured
in wet labs, they exist for tens of thousands of compound/isoform pairs,
and they do not depend on anyone's structural interpretation.

If the mined variable positions genuinely determine isoform selectivity,
then a specific prediction follows:

  isoform pairs differing at MORE mined contact positions should show
  LARGER experimental selectivity spreads across compounds measured
  against both.

That is one pre-specified statistic computed over all isoform pairs at
once. It cannot be cherry-picked a residue at a time.

CONTROLS
--------
A positive correlation alone is not enough -- isoform pairs that diverge
at contact positions also diverge everywhere else, so the result could
merely say "different proteins bind differently." Two controls separate
the hypotheses:

  C1  whole-protein sequence divergence. If overall divergence predicts
      selectivity just as well, the mined positions add nothing.
  C2  random size-matched sets of NON-contact positions, resampled to
      build a null distribution. If the real contact set does not beat
      that null, the tool is finding nothing special.

Under H2, the mined positions should perform no better than C1/C2.
Under H1, they should beat both.
"""

import itertools
import json
import os
import random
import sys

import numpy as np
import requests
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from posegate.selectivity import fetch_uniprot_sequence, build_alignment_map, top_residue_numbers

CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
CACHE_DIR = "data/chembl_cache"
REFERENCE = 'CA2'
TOP_N = 10
N_NULL_DRAWS = 2000
RANDOM_SEED = 0

# ChEMBL target IDs resolved from UniProt accession via the ChEMBL target
# API, not guessed -- an initial guess of CHEMBL3510 for CA13 turned out
# to be Carbonic anhydrase 14.
ISOFORMS = {
    'CA1':  {'acc': 'P00915', 'chembl': 'CHEMBL261',  'mined': 'data/ca1_verified/mined_result.json'},
    'CA2':  {'acc': 'P00918', 'chembl': 'CHEMBL205',  'mined': 'data/ca_verified/mined_result.json'},
    'CA4':  {'acc': 'P22748', 'chembl': 'CHEMBL3729', 'mined': 'data/ca4_verified/mined_result.json'},
    'CA7':  {'acc': 'P43166', 'chembl': 'CHEMBL2326', 'mined': 'data/ca7_verified/mined_result.json'},
    'CA9':  {'acc': 'Q16790', 'chembl': 'CHEMBL3594', 'mined': 'data/ca9_verified/mined_result.json'},
    'CA12': {'acc': 'O43570', 'chembl': 'CHEMBL3242', 'mined': 'data/ca12_verified/mined_result.json'},
    'CA13': {'acc': 'Q8N1Q1', 'chembl': 'CHEMBL3912', 'mined': 'data/ca13_verified/mined_result.json'},
}
MIN_SHARED_COMPOUNDS = 25


def fetch_ki(name, target_id):
    """{molecule_chembl_id: median pChEMBL} for one target, cached.
    pChEMBL is ChEMBL's own normalised -log10(molar) value, so Ki, IC50
    and Kd are already on one comparable scale and unit-parsing errors
    are not this script's problem. A compound measured several times is
    reduced to its median rather than an arbitrary single record."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = os.path.join(CACHE_DIR, f"{name}.json")
    if os.path.exists(path):
        with open(path) as f:
            return {k: v for k, v in json.load(f).items()}

    raw = {}
    offset, limit = 0, 1000
    while True:
        r = requests.get(CHEMBL_ACTIVITY, timeout=120, params={
            'target_chembl_id': target_id, 'standard_type': 'Ki',
            'pchembl_value__isnull': 'false', 'limit': limit, 'offset': offset})
        r.raise_for_status()
        page = r.json()
        for a in page['activities']:
            mol, val = a.get('molecule_chembl_id'), a.get('pchembl_value')
            if mol and val is not None:
                raw.setdefault(mol, []).append(float(val))
        if not page['page_meta'].get('next'):
            break
        offset += limit
    merged = {m: float(np.median(v)) for m, v in raw.items()}
    with open(path, 'w') as f:
        json.dump(merged, f)
    return merged


def mined_contact_positions(sequences, maps_to_ref):
    """Reference-coordinate positions that are a top-N conserved contact
    in at least one isoform -- the tool's actual output, unchanged."""
    positions = set()
    for name, info in ISOFORMS.items():
        for native in top_residue_numbers(info['mined'], TOP_N):
            if native > len(sequences[name]):
                continue  # retained metal ion, no sequence position
            ref_pos = native if name == REFERENCE else maps_to_ref[name].get(native)
            if ref_pos is not None:
                positions.add(ref_pos)
    return sorted(positions)


def divergence(a, b, positions, sequences, maps_from_ref):
    """How many of `positions` (reference coordinates) carry different
    amino acids in isoforms a and b. Positions unaligned in either are
    skipped, not counted as differences."""
    n = 0
    for p in positions:
        aa_a = sequences[a][p - 1] if a == REFERENCE else (
            sequences[a][maps_from_ref[a][p] - 1] if p in maps_from_ref[a] else None)
        aa_b = sequences[b][p - 1] if b == REFERENCE else (
            sequences[b][maps_from_ref[b][p] - 1] if p in maps_from_ref[b] else None)
        if aa_a and aa_b and aa_a != aa_b:
            n += 1
    return n


def main():
    random.seed(RANDOM_SEED)
    names = sorted(ISOFORMS)

    print("Fetching experimental Ki data from ChEMBL (cached after first run)...")
    ki = {}
    for name in names:
        ki[name] = fetch_ki(name, ISOFORMS[name]['chembl'])
        print(f"  {name}: {len(ki[name])} compounds with pChEMBL")

    sequences = {n: fetch_uniprot_sequence(ISOFORMS[n]['acc']) for n in names}
    ref_seq = sequences[REFERENCE]
    maps_from_ref, maps_to_ref = {}, {}
    for n in names:
        if n == REFERENCE:
            continue
        m = build_alignment_map(ref_seq, sequences[n])
        maps_from_ref[n] = m
        maps_to_ref[n] = {v: k for k, v in m.items()}

    contact_positions = mined_contact_positions(sequences, maps_to_ref)
    non_contact = [p for p in range(1, len(ref_seq) + 1) if p not in set(contact_positions)]
    print(f"\nMined contact positions (the tool's output): {len(contact_positions)}")
    print(f"  {contact_positions}")

    rows = []
    for a, b in itertools.combinations(names, 2):
        shared = set(ki[a]) & set(ki[b])
        if len(shared) < MIN_SHARED_COMPOUNDS:
            print(f"  skipping {a}/{b}: only {len(shared)} shared compounds")
            continue
        deltas = [abs(ki[a][m] - ki[b][m]) for m in shared]
        rows.append({
            'pair': f"{a}/{b}", 'a': a, 'b': b, 'n_shared': len(shared),
            'observed': float(np.median(deltas)),
            'predicted': divergence(a, b, contact_positions, sequences, maps_from_ref),
            'whole_protein': divergence(a, b, list(range(1, len(ref_seq) + 1)),
                                         sequences, maps_from_ref),
        })

    print(f"\n{'pair':<12}{'n_cmpd':>8}{'contact div':>13}{'whole-prot div':>16}{'median |dpKi|':>15}")
    print('-' * 64)
    for r in sorted(rows, key=lambda r: -r['predicted']):
        print(f"{r['pair']:<12}{r['n_shared']:>8}{r['predicted']:>13}"
              f"{r['whole_protein']:>16}{r['observed']:>15.3f}")

    pred = [r['predicted'] for r in rows]
    obs = [r['observed'] for r in rows]
    whole = [r['whole_protein'] for r in rows]

    rho, p = spearmanr(pred, obs)
    rho_c1, p_c1 = spearmanr(whole, obs)

    print(f"\n{'='*64}")
    print(f"MAIN TEST   mined contact divergence vs experimental selectivity")
    print(f"            Spearman rho = {rho:+.3f}   p = {p:.4f}   (n = {len(rows)} pairs)")
    print(f"\nCONTROL C1  whole-protein divergence vs experimental selectivity")
    print(f"            Spearman rho = {rho_c1:+.3f}   p = {p_c1:.4f}")

    # C2: null distribution from random size-matched non-contact sets.
    k = len(contact_positions)
    null = []
    for _ in range(N_NULL_DRAWS):
        sample = random.sample(non_contact, k)
        d = [divergence(r['a'], r['b'], sample, sequences, maps_from_ref) for r in rows]
        if len(set(d)) > 1:
            null.append(spearmanr(d, obs)[0])
    null = np.array([x for x in null if not np.isnan(x)])
    better = float(np.mean(null >= rho)) if len(null) else float('nan')

    print(f"\nCONTROL C2  {len(null)} random size-matched NON-contact position sets")
    print(f"            null rho: mean {null.mean():+.3f}, sd {null.std():.3f}, "
          f"95th pct {np.percentile(null, 95):+.3f}")
    print(f"            fraction of random sets matching or beating the real one: "
          f"{better:.4f}")
    print(f"{'='*64}")

    with open('data/selectivity_vs_experiment.json', 'w') as f:
        json.dump({'rows': rows, 'contact_positions': contact_positions,
                    'rho': rho, 'p': p, 'rho_whole_protein': rho_c1,
                    'p_whole_protein': p_c1, 'null_empirical_p': better,
                    'n_null_draws': int(len(null))}, f, indent=2)
    print("\nWrote data/selectivity_vs_experiment.json")


if __name__ == "__main__":
    main()
