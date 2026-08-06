# posegate/tests/test_prep_ensemble.py
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from prep_ensemble import filter_chains


def test_filter_chains_drops_ter_records_for_filtered_out_chains(tmp_path):
    """A dropped chain's TER record must be dropped along with its
    ATOM/HETATM lines, not passed through -- an unmatched TER (no
    preceding ATOM in that model) crashes OpenMM's PDB parser outright.
    Found on chymotrypsin, a 3-chain deposition where an earlier version
    of this function wrote every filtered-out chain's TER through
    unconditionally."""
    pdb_in = tmp_path / "in.pdb"
    pdb_out = tmp_path / "out.pdb"
    pdb_in.write_text(
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "TER       2      ALA A   1\n"
        "ATOM      3  N   GLY B   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "TER       4      GLY B   1\n"
        "END\n"
    )

    filter_chains(str(pdb_in), keep_chains={'A'}, out_path=str(pdb_out))

    lines = pdb_out.read_text().splitlines()
    assert any(l.startswith('ATOM') and l[21] == 'A' for l in lines)
    assert not any(l.startswith('ATOM') and l[21] == 'B' for l in lines)
    assert not any(l.startswith('TER') and l[21] == 'B' for l in lines)
    assert any(l.startswith('TER') and l[21] == 'A' for l in lines)
    assert any(l.startswith('END') for l in lines)  # non-chain-scoped lines pass through


def test_filter_chains_keeps_multiple_requested_chains(tmp_path):
    pdb_in = tmp_path / "in.pdb"
    pdb_out = tmp_path / "out.pdb"
    pdb_in.write_text(
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "TER       2      ALA A   1\n"
        "ATOM      3  N   GLY B   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "TER       4      GLY B   1\n"
        "ATOM      5  N   SER C   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "TER       6      SER C   1\n"
    )

    filter_chains(str(pdb_in), keep_chains={'A', 'C'}, out_path=str(pdb_out))

    lines = pdb_out.read_text().splitlines()
    chains_present = {l[21] for l in lines if l.startswith(('ATOM', 'TER'))}
    assert chains_present == {'A', 'C'}
