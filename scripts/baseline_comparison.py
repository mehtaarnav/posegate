# posegate/scripts/baseline_comparison.py
"""Does the sophisticated method (specific-interaction filtering + LOO
self-validation) actually change which residues rank at the top,
compared to the simplest possible baseline: rank by raw frequency of
ANY contact (including the unfiltered VdWContact superset) across the
ensemble, no self-validation at all?

Run against the two already-confirmed cases (CA2/CA9's Phe131/Val,
CDK2/CDK9's Leu83/Cys106 -- see conversation) using data already mined,
no new fetching or compute. If the naive baseline produces the same
top-N, the sophistication isn't earning its keep for THESE SPECIFIC
discoveries (though it may still matter for confidence quantification,
which a naive frequency count has no equivalent of). If it produces a
different, worse top-N, that's real evidence the filtering/validation
is doing something, not just adding complexity.
"""

import json


def naive_top_n(mined_path: str, top_n: int) -> list:
    """Baseline: rank by raw frequency across ALL interaction rows,
    including VdWContact (which the sophisticated method deliberately
    excludes as an uninformative superset -- see conserved_contacts.py's
    module docstring), with no LOO self-validation step at all."""
    with open(mined_path) as f:
        data = json.load(f)
    best_freq = {}
    for row in data['mined']:
        num = int(''.join(c for c in row['residue'] if c.isdigit()))
        best_freq[num] = max(best_freq.get(num, 0), row['frequency'])
    ranked = sorted(best_freq.items(), key=lambda kv: (-kv[1], kv[0]))
    return [num for num, _ in ranked[:top_n]]


def sophisticated_top_n(mined_path: str, top_n: int) -> list:
    """The method actually used elsewhere in this project: VdWContact
    excluded, ranked by frequency of specific interaction types only."""
    with open(mined_path) as f:
        data = json.load(f)
    seen = []
    for row in data['mined']:
        if row['interaction'] == 'VdWContact':
            continue
        num = int(''.join(c for c in row['residue'] if c.isdigit()))
        if num not in seen:
            seen.append(num)
        if len(seen) >= top_n:
            break
    return seen


CASES = [
    {
        'name': 'CA2 vs CA9 (confirmed: Phe130 CA2 / Val262 CA9, raw UniProt numbering, '
                'each isoform in its own native numbering -- see ca_family_selectivity.py)',
        'targets': {'CA2': 130, 'CA9': 262},
        'files': {'CA2': 'data/ca_verified/mined_result.json',
                  'CA9': 'data/ca9_verified/mined_result.json'},
    },
    {
        'name': 'CDK2 vs CDK9 (confirmed: Leu83 CDK2 / Cys106 CDK9, each in native numbering)',
        'targets': {'CDK2': 83, 'CDK9': 106},
        'files': {'CDK2': 'data/cdk2_verified/mined_result.json',
                  'CDK9': 'data/cdk9_verified/mined_result.json'},
    },
]


def main():
    top_n = 10
    for case in CASES:
        print(f"\n{'='*70}\n{case['name']}\n{'='*70}")
        for isoform, path in case['files'].items():
            naive = naive_top_n(path, top_n)
            sophisticated = sophisticated_top_n(path, top_n)
            target = case['targets'][isoform]
            print(f"\n{isoform}:")
            print(f"  naive (raw freq incl. VdW, no LOO):        {naive}")
            print(f"  sophisticated (specific interactions, LOO): {sophisticated}")
            print(f"  target position {target} in naive top-{top_n}?        "
                  f"{'YES at rank ' + str(naive.index(target)+1) if target in naive else 'NO'}")
            print(f"  target position {target} in sophisticated top-{top_n}? "
                  f"{'YES at rank ' + str(sophisticated.index(target)+1) if target in sophisticated else 'NO'}")
            print(f"  same top-{top_n} set: {set(naive) == set(sophisticated)}")


if __name__ == "__main__":
    main()
