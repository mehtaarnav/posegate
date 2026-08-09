# posegate/scripts/ligand_conditioned_test.py
"""Ligand-conditioned test of the selectivity claim.

The pair-level test (selectivity_vs_experiment.py) failed: rho = +0.244,
p = 0.287. The variance decomposition then showed why -- 74.5% of dpKi
variance sits WITHIN isoform pairs, driven by which compound is used,
and only 25.5% BETWEEN pairs. Collapsing each pair to a median discarded
three quarters of the signal by construction.

This conditions on the ligand instead of averaging over it. For a
compound whose crystal structure is in the mined ensemble, we know
exactly which residue positions THAT compound touches. So:

  PRE-SPECIFIED HYPOTHESIS
  For a compound c with a structure in isoform A, and any isoform B
  where c also has measured Ki, the number of positions c contacts in A
  that DIFFER between A and B should predict |pKi(c,A) - pKi(c,B)|.

  A compound that touches nothing which differs between two isoforms has
  no structural basis to discriminate them; one that touches several
  divergent positions does.

  PRE-SPECIFIED CONTROL
  The number of positions c contacts that are IDENTICAL between A and B
  should NOT predict selectivity. If it does equally well, the signal is
  "compounds with more contacts overall are more selective" -- a size
  artifact, not a selectivity mechanism.

Both statistics are fixed before the run. Both are reported regardless
of outcome.
"""

import glob
import json
import os
import sys

import numpy as np
import requests
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from posegate.selectivity import fetch_uniprot_sequence, build_alignment_map
from posegate.conserved_contacts import _structure_contact_residues
from selectivity_vs_experiment import ISOFORMS, fetch_ki, REFERENCE

UNICHEM = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
CACHE = "data/ligand_conditioned_cache.json"


def pdb_to_chembl(codes):
    """{PDB chem-comp id: ChEMBL id} via UniChem. src id 3 is the PDB,
    id 1 is ChEMBL. Cached -- this is ~70 network calls otherwise."""
    cache_path = "data/pdb_to_chembl.json"
    cached = {}
    if os.path.exists(cache_path):
        cached = json.load(open(cache_path))
    out = {c: cached[c] for c in codes if c in cached and cached[c]}
    todo = [c for c in codes if c not in out]
    for c in todo:
        chembl = None
        try:
            r = requests.post(UNICHEM, timeout=30,
                              json={'type': 'sourceID', 'compound': c, 'sourceID': 3})
            if r.status_code == 200:
                for comp in r.json().get('compounds', []):
                    for src in comp.get('sources', []):
                        if src.get('id') == 1:
                            chembl = src.get('compoundId')
                            break
                    if chembl:
                        break
        except Exception:
            pass
        cached[c] = chembl
        if chembl:
            out[c] = chembl
    json.dump(cached, open(cache_path, 'w'))
    return out


def structure_ligand_code(prepped_dir, pdb_id):
    path = os.path.join(prepped_dir, f"{pdb_id}_ligand_raw.pdb")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            if line.startswith('HETATM'):
                return line[17:20].strip()
    return None


def main():
    names = sorted(ISOFORMS)
    sequences = {n: fetch_uniprot_sequence(ISOFORMS[n]['acc']) for n in names}
    ref_seq = sequences[REFERENCE]
    from_ref, to_ref = {}, {}
    for n in names:
        if n == REFERENCE:
            continue
        m = build_alignment_map(ref_seq, sequences[n])
        from_ref[n] = m
        to_ref[n] = {v: k for k, v in m.items()}

    def aa_at(iso, ref_pos):
        if iso == REFERENCE:
            return ref_seq[ref_pos - 1] if ref_pos <= len(ref_seq) else None
        p = from_ref[iso].get(ref_pos)
        return sequences[iso][p - 1] if p else None

    # Per-structure contacts, cached (recomputing IFPs is the slow part).
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    structures = []
    for iso in names:
        d = os.path.dirname(ISOFORMS[iso]['mined'])
        for pkl in sorted(glob.glob(os.path.join(d, '*_receptor_h.pkl'))):
            pdb_id = os.path.basename(pkl).replace('_receptor_h.pkl', '')
            sdf = os.path.join(d, f'{pdb_id}_ligand_h.sdf')
            if not os.path.exists(sdf):
                continue
            key = f"{iso}:{pdb_id}"
            if key not in cache:
                res = _structure_contact_residues(
                    {'pdb_id': pdb_id, 'ligand_sdf': sdf, 'receptor_pdb': pkl})
                positions = []
                if res:
                    for r in res:
                        digits = ''.join(ch for ch in r if ch.isdigit())
                        if not digits:
                            continue
                        native = int(digits)
                        rp = native if iso == REFERENCE else to_ref[iso].get(native)
                        if rp:
                            positions.append(rp)
                cache[key] = sorted(set(positions))
            lig = structure_ligand_code(d, pdb_id)
            if lig and cache[key]:
                structures.append({'iso': iso, 'pdb_id': pdb_id, 'ligand': lig,
                                    'positions': cache[key]})
    json.dump(cache, open(CACHE, 'w'))
    print(f"Structures with usable contacts: {len(structures)}")

    codes = sorted({s['ligand'] for s in structures})
    mapping = pdb_to_chembl(codes)
    print(f"Ligand codes mapped to ChEMBL: {len(mapping)}/{len(codes)}")

    ki = {n: fetch_ki(n, ISOFORMS[n]['chembl']) for n in names}

    rows = []
    for s in structures:
        chembl = mapping.get(s['ligand'])
        if not chembl:
            continue
        A = s['iso']
        if chembl not in ki[A]:
            continue
        for B in names:
            if B == A or chembl not in ki[B]:
                continue
            divergent = invariant = 0
            for p in s['positions']:
                a, b = aa_at(A, p), aa_at(B, p)
                if a is None or b is None:
                    continue
                if a != b:
                    divergent += 1
                else:
                    invariant += 1
            rows.append({'compound': chembl, 'ligand': s['ligand'], 'pdb_id': s['pdb_id'],
                          'A': A, 'B': B, 'n_contacts': len(s['positions']),
                          'divergent': divergent, 'invariant': invariant,
                          'observed': abs(ki[A][chembl] - ki[B][chembl])})

    print(f"Usable (compound, isoform-pair) observations: {len(rows)}")
    if len(rows) < 20:
        print("Too few observations for a meaningful correlation -- reporting and stopping.")
        for r in rows:
            print(' ', r)
        return

    div = [r['divergent'] for r in rows]
    inv = [r['invariant'] for r in rows]
    obs = [r['observed'] for r in rows]

    rho_d, p_d = spearmanr(div, obs)
    rho_i, p_i = spearmanr(inv, obs)

    print(f"\n{'='*68}")
    print(f"MAIN     divergent contacts vs |dpKi|   rho = {rho_d:+.3f}  p = {p_d:.4f}")
    print(f"CONTROL  invariant contacts vs |dpKi|   rho = {rho_i:+.3f}  p = {p_i:.4f}")
    print(f"         n = {len(rows)} observations, "
          f"{len({r['compound'] for r in rows})} distinct compounds")
    print(f"{'='*68}")

    print(f"\nDivergent-contact count distribution: "
          f"{sorted(set(div))} (n distinct = {len(set(div))})")
    print(f"Median |dpKi| by divergent-contact count:")
    for d in sorted(set(div)):
        vals = [r['observed'] for r in rows if r['divergent'] == d]
        print(f"  {d} divergent contacts: n={len(vals):<5} median |dpKi| = {np.median(vals):.2f}")

    json.dump({'rows': rows, 'rho_divergent': rho_d, 'p_divergent': p_d,
                'rho_invariant': rho_i, 'p_invariant': p_i},
              open('data/ligand_conditioned_result.json', 'w'), indent=2)
    print("\nWrote data/ligand_conditioned_result.json")


if __name__ == "__main__":
    main()
