# posegate/scripts/ca_family_selectivity.py
"""Cross-paralog selectivity comparison for the human carbonic anhydrase
family: CA II (P00918, the validated reference isoform), CA IX (Q16790,
tumor-associated), CA XII (O43570, tumor-associated). Logic lives in
posegate.selectivity; this script is just the CA-specific configuration
and report formatting.

CA II is the well-characterized, non-tumor "off-target" that most CA
inhibitors bind nonselectively; CA IX/XII selectivity over CA II is a
real, long-standing drug design problem.

Verification status of this run's specific top-N output (checked against
independent literature, not just self-consistency -- see conversation):
  - CA2-ref 130 (raw UniProt numbering; classical "CA II numbering" 131):
    CONFIRMED. Phe in CA2, Val in CA9 -- matches PMC7534198's documented
    F131V CA IX-selectivity mechanism exactly, including the classical-
    numbering +1 offset from raw UniProt indexing used throughout this
    module.
  - CA2-ref 131 (classical 132): CONFIRMED. Gly in CA2, Asp in CA9 --
    matches published hCA II/hCA IX active-site rim literature ("Gly in
    hCA II and Asp in hCA IX" at this position) exactly.
  - CA2-ref 200 (classical 201, CA12-unique candidate): NOT confirmed.
    Literature places Pro201/Pro202 as part of the broadly conserved
    hydrophilic rim, present in CA9's active site too -- not a CA12-
    specific feature. Consistent with this module's own data: the
    residue is Pro in both CA2 and CA12 (no chemical divergence at all),
    so its appearance in the "unique to CA12" bucket here is a frequency-
    ranking artifact, not a real selectivity signal.
"""

from posegate.selectivity import compare_isoforms

ISOFORMS = {
    'CA2':  {'acc': 'P00918', 'mined': 'data/ca_verified/mined_result.json'},
    'CA9':  {'acc': 'Q16790', 'mined': 'data/ca9_verified/mined_result.json'},
    'CA12': {'acc': 'O43570', 'mined': 'data/ca12_verified/mined_result.json'},
}
REFERENCE = 'CA2'
TOP_N = 10


def main():
    result = compare_isoforms(ISOFORMS, REFERENCE, TOP_N)

    for name, length in result['sequences'].items():
        print(f"{name} ({ISOFORMS[name]['acc']}): {length} residues")

    print(f"\nTop-{TOP_N} mined residues (own numbering, own protein):")
    for name, nums in result['top_native'].items():
        print(f"  {name}: {nums}")

    for name, missed in result['unmapped'].items():
        print(f"\n{name}: {len(missed)} top residue(s) fell in a gapped/unaligned region "
              f"and have no {REFERENCE}-coordinate equivalent: {missed}")

    print(f"\nTop-{TOP_N} residues, all expressed in {REFERENCE} reference numbering:")
    for name, positions in result['top_in_reference_coords'].items():
        print(f"  {name}: {sorted(positions)}")

    print(f"\nShared across all isoforms (pan-CA catalytic scaffold candidates): "
          f"{sorted(result['shared_by_all'])}")
    for (a, b), positions in result['pairwise_only'].items():
        others = [n for n in ISOFORMS if n not in (a, b)]
        print(f"Shared by {a}+{b} only (not {'+'.join(others)}): {sorted(positions)}")
    for name, positions in result['unique_to'].items():
        print(f"Unique to {name} (selectivity-relevant candidates): {sorted(positions)}")


if __name__ == "__main__":
    main()
