# posegate/scripts/compare_miners.py
"""Quantitative comparison of posegate's ProLIF-based conserved-contact
miner against the same mining done from PLIP's output, across the 5
protein families. Interaction-type vocabularies differ between the two
tools (ProLIF splits HBDonor/HBAcceptor by ligand role and distinguishes
FaceToFace/EdgeToFace pi-stacking; PLIP reports one HBond bucket and one
PiStack bucket) so comparison is done at the residue level: does each
tool flag this residue as *any* kind of specific (non-VdW) contact, and
at what frequency."""

import argparse
import json
import re

# ProLIF's VdWContact is the loosest possible criterion (any atoms within
# vdW radii) and has no PLIP equivalent (PLIP only reports interactions it
# considers specific), so it's excluded from the residue-level comparison
# to avoid making posegate look more "complete" by a criterion PLIP was
# never trying to report in the first place.
POSEGATE_SPECIFIC = {'HBDonor', 'HBAcceptor', 'Hydrophobic', 'FaceToFace', 'EdgeToFace'}


def normalize_residue(residue: str, source: str) -> str:
    """Normalize to 'number.chain'. posegate residues come from ProLIF's
    ResidueId string repr, e.g. 'ASN140.A' (name+number+chain); PLIP's are
    already 'number.chain' (its resnr/reschain fields have no residue-name
    prefix)."""
    if source == 'posegate':
        m = re.match(r'^[A-Za-z]*(\d+)\.(.+)$', residue)
        if m:
            return f"{m.group(1)}.{m.group(2)}"
    return residue


def residue_set_chain_aware(freq_table: list, source: str) -> dict:
    out = {}
    for r in freq_table:
        if source == 'posegate' and r['interaction'] not in POSEGATE_SPECIFIC:
            continue
        residue = normalize_residue(r['residue'], source)
        out[residue] = max(out.get(residue, 0.0), r['frequency'])
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--families", nargs='+', required=True,
                         help="Family names matching plip_out/{name}_posegate_freq.json and _plip_freq.json")
    parser.add_argument("--freq_threshold", type=float, default=0.5,
                         help="Minimum frequency to count a residue as 'conserved'")
    args = parser.parse_args()

    print(f"{'Family':<12}{'posegate conserved':>20}{'PLIP conserved':>18}{'agree':>8}{'posegate-only':>15}{'PLIP-only':>12}")
    print('-' * 90)

    for family in args.families:
        with open(f"plip_out/{family}_posegate_freq.json") as f:
            pg_table = json.load(f)
        with open(f"plip_out/{family}_plip_freq.json") as f:
            plip_table = json.load(f)

        pg_res = residue_set_chain_aware(pg_table, 'posegate')
        plip_res = residue_set_chain_aware(plip_table, 'plip')

        pg_conserved = {r for r, f in pg_res.items() if f >= args.freq_threshold}
        plip_conserved = {r for r, f in plip_res.items() if f >= args.freq_threshold}

        agree = pg_conserved & plip_conserved
        pg_only = pg_conserved - plip_conserved
        plip_only = plip_conserved - pg_conserved

        print(f"{family:<12}{len(pg_conserved):>20}{len(plip_conserved):>18}{len(agree):>8}{len(pg_only):>15}{len(plip_only):>12}")
        print(f"  posegate conserved: {sorted(pg_conserved)}")
        print(f"  PLIP conserved:     {sorted(plip_conserved)}")
        print(f"  agree:              {sorted(agree)}")
        print()


if __name__ == "__main__":
    main()
