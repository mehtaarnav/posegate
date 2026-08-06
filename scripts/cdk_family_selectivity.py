# posegate/scripts/cdk_family_selectivity.py
"""Cross-paralog selectivity comparison for CDK2 (P24941) vs CDK9
(P50750) -- the second family this method was run on, to check whether
the CA result (see ca_family_selectivity.py) was a real capability or a
fluke on one enzyme class. Logic lives in posegate.selectivity; this
script is just the CDK-specific configuration and report formatting.

CDK2 vs CDK9 selectivity is a real, actively pursued oncology problem:
CDK9 (transcriptional, via P-TEFb) inhibitors are being developed for
cancer, and off-target CDK2 (cell-cycle) inhibition is a real selectivity
liability researchers explicitly try to avoid.

Verification status (checked against independent literature -- see
conversation): CDK2's top-1 mined residue (Leu83) and CDK9's top-1 mined
residue (Cys106) align to the exact same physical hinge position.
CONFIRMED verbatim against literature: "The hinge region residue Cys106
in CDK9 (corresponding to Leu83 in CDK2)... is a key strategy for
achieving selective inhibition of CDK9 over CDK2" -- exact residue
numbers, exact correspondence direction, no numbering-convention
correction needed (unlike the CA case).
"""

from posegate.selectivity import compare_isoforms

ISOFORMS = {
    'CDK2': {'acc': 'P24941', 'mined': 'data/cdk2_verified/mined_result.json'},
    'CDK9': {'acc': 'P50750', 'mined': 'data/cdk9_verified/mined_result.json'},
}
REFERENCE = 'CDK2'
TOP_N = 10


def main():
    result = compare_isoforms(ISOFORMS, REFERENCE, TOP_N)

    for name, length in result['sequences'].items():
        print(f"{name} ({ISOFORMS[name]['acc']}): {length} residues")

    print(f"\nTop-{TOP_N} mined residues (own numbering, own protein):")
    for name, nums in result['top_native'].items():
        print(f"  {name}: {nums}")

    for name, missed in result['unmapped'].items():
        print(f"\n{name}: {len(missed)} top residue(s) unmapped (gapped region): {missed}")

    print(f"\nTop-{TOP_N} residues, all expressed in {REFERENCE} reference numbering:")
    for name, positions in result['top_in_reference_coords'].items():
        print(f"  {name}: {sorted(positions)}")

    print(f"\nShared (conserved kinase-domain scaffold candidates): {sorted(result['shared_by_all'])}")
    for name, positions in result['unique_to'].items():
        print(f"Unique to {name} (selectivity-relevant candidates): {sorted(positions)}")


if __name__ == "__main__":
    main()
