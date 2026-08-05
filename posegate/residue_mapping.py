# posegate/posegate/residue_mapping.py
"""Cross-structure residue-number remapping via PDBe SIFTS.

conserved_contacts.py keys every mined contact on (chain, author residue
number) as read directly off each structure's PDB file, on the implicit
assumption that the same (chain, resnum) pair means the same physical
residue across every structure in an ensemble. That assumption breaks
whenever an ensemble is auto-fetched (scripts/mine_target.py) rather than
hand-curated: different depositions of the same protein can use different
author numbering (different construct boundaries, tags, isoforms). Found
directly on ERalpha: position 353 is GLU in one fetched structure and LYS
in another -- not a numbering offset, a different residue entirely at the
same author number -- which fragmented the conserved-contact vote across
several numbering registers and drove leave-one-out top-1 accuracy to 0%.

This remaps every structure's author residue numbers onto UniProt residue
numbers before mining, using PDBe's own SIFTS PDB-to-UniProt mapping
(https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/), so (chain, resnum)
means the same physical position across every structure regardless of
each deposition's own author numbering.
"""

from typing import Dict, List, Tuple

import requests

SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
RESIDUE_LISTING_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id}/chain/{chain}"


def fetch_uniprot_segments(pdb_id: str, uniprot_acc: str) -> List[Dict]:
    """SIFTS's PDB-to-UniProt segments for one accession in one PDB entry,
    in PDBe's own sequential 'residue_number' coordinate (SEQRES-based,
    always present and monotonic) rather than author numbering:
    [{'chain', 'label_start', 'label_end', 'unp_start'}, ...].

    Deliberately not author_residue_number: SIFTS reports that as null at
    segment boundaries whenever a segment's author numbering isn't
    provably linear (seen directly on 3ERT), so a caller assuming a fixed
    author-number offset per segment can silently mismap residues past
    any internal insertion/renumbering. residue_number has no such gap;
    the corresponding author number for each residue_number is looked up
    per-structure instead, in build_residue_map, via residue_listing.

    Raises ValueError if this entry's SIFTS mapping doesn't include that
    accession at all (wrong accession, or an entry SIFTS hasn't mapped)."""
    resp = requests.get(SIFTS_URL.format(pdb_id=pdb_id.lower()), timeout=30)
    resp.raise_for_status()
    data = resp.json()
    entry = data.get(pdb_id.lower(), {})
    uniprot_block = entry.get('UniProt', {})
    if uniprot_acc not in uniprot_block:
        raise ValueError(f"{pdb_id}: SIFTS has no mapping to UniProt {uniprot_acc} "
                          f"(mapped accessions here: {list(uniprot_block.keys())})")
    segments = []
    for m in uniprot_block[uniprot_acc]['mappings']:
        segments.append({
            'chain': m['chain_id'],
            'label_start': m['start']['residue_number'],
            'label_end': m['end']['residue_number'],
            'unp_start': m['unp_start'],
        })
    return segments


def fetch_author_numbers(pdb_id: str, chain: str) -> Dict[int, int]:
    """{residue_number (PDBe sequential label): author_residue_number}
    for one chain of one structure, straight from PDBe's residue listing
    -- the actual numbers written in the deposited PDB file's ATOM
    records, which is what conserved_contacts.py reads."""
    resp = requests.get(RESIDUE_LISTING_URL.format(pdb_id=pdb_id.lower(), chain=chain), timeout=30)
    resp.raise_for_status()
    data = resp.json()
    residues = data[pdb_id.lower()]['molecules'][0]['chains'][0]['residues']
    return {r['residue_number']: r['author_residue_number'] for r in residues
            if r['author_residue_number'] is not None}


def build_residue_map(pdb_id: str, uniprot_acc: str) -> Dict[Tuple[str, int], int]:
    """{(chain, author_resnum): uniprot_resnum} for one structure, built
    per-residue (not by assuming a linear author-number offset across a
    whole segment -- see fetch_uniprot_segments for why that broke on
    3ERT). A residue_number covered by a UniProt segment but missing from
    the residue listing (usually because it was never observed in the
    density) is simply absent from the returned map, not mismapped."""
    segments = fetch_uniprot_segments(pdb_id, uniprot_acc)
    mapping: Dict[Tuple[str, int], int] = {}
    author_by_chain: Dict[str, Dict[int, int]] = {}
    for seg in segments:
        chain = seg['chain']
        if chain not in author_by_chain:
            author_by_chain[chain] = fetch_author_numbers(pdb_id, chain)
        auth_lookup = author_by_chain[chain]
        for label_num in range(seg['label_start'], seg['label_end'] + 1):
            auth_num = auth_lookup.get(label_num)
            if auth_num is None:
                continue
            mapping[(chain, auth_num)] = seg['unp_start'] + (label_num - seg['label_start'])
    return mapping


def remap_residue_numbers(pdb_path: str, out_path: str,
                           residue_map: Dict[Tuple[str, int], int]) -> int:
    """Rewrites ATOM/HETATM resSeq fields (PDB columns 23-26, 1-indexed)
    to UniProt residue numbers wherever (chain, author_resnum) is in
    residue_map. Lines outside the map (ligand HETATM, unmapped chains,
    residues SIFTS didn't cover) pass through with their original number
    unchanged, not dropped -- an unmapped receptor residue still needs to
    exist for van der Waals contact detection, it's just not comparable
    across structures by number.

    Returns the count of atom lines actually remapped, so a caller can
    detect a near-total miss (e.g. a wrong chain ID convention) instead
    of silently mining an ensemble that's still inconsistent.
    """
    n_remapped = 0
    with open(pdb_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM')):
                chain = line[21]
                try:
                    auth_num = int(line[22:26])
                except ValueError:
                    fout.write(line)
                    continue
                key = (chain, auth_num)
                if key in residue_map:
                    line = line[:22] + f"{residue_map[key]:>4}" + line[26:]
                    n_remapped += 1
            fout.write(line)
    return n_remapped
