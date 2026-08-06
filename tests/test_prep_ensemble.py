# posegate/tests/test_prep_ensemble.py
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from prep_ensemble import filter_chains, chain_residue_sequences, is_asymmetric_multichain


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


def _residue_line(chain, resnum, resname, serial=1):
    return (f"ATOM  {serial:>5}  CA  {resname} {chain}{resnum:>4}      "
            f"0.000   0.000   0.000  1.00  0.00           C\n")


def test_is_asymmetric_multichain_false_for_a_true_symmetric_homodimer(tmp_path):
    """HIV protease's homodimer: both chains identical sequence -- any
    cross-chain residue-label mixing is harmless, must NOT be flagged."""
    pdb = tmp_path / "homodimer.pdb"
    lines = []
    for i, resname in enumerate(['ALA', 'GLY', 'SER'], start=1):
        lines.append(_residue_line('A', i, resname))
        lines.append(_residue_line('B', i, resname))
    pdb.write_text(''.join(lines))

    assert is_asymmetric_multichain(str(pdb), keep_chains={'A', 'B'}) is False


def test_is_asymmetric_multichain_true_for_non_identical_chains(tmp_path):
    """Chymotrypsin-like: three non-identical fragments of one cleaved
    polypeptide -- must be flagged, this is the real risk case."""
    pdb = tmp_path / "chymotrypsin_like.pdb"
    lines = [
        _residue_line('A', 1, 'ILE'), _residue_line('A', 2, 'VAL'),
        _residue_line('B', 1, 'GLY'), _residue_line('B', 2, 'GLY'),
        _residue_line('C', 1, 'SER'), _residue_line('C', 2, 'ASP'),
    ]
    pdb.write_text(''.join(lines))

    assert is_asymmetric_multichain(str(pdb), keep_chains={'A', 'B', 'C'}) is True


def test_is_asymmetric_multichain_false_for_a_single_chain(tmp_path):
    """CA2/CDK2-like: one chain near the ligand -- nothing to be
    inconsistent with, must NOT be flagged."""
    pdb = tmp_path / "single_chain.pdb"
    pdb.write_text(_residue_line('A', 1, 'ALA') + _residue_line('A', 2, 'GLY'))

    assert is_asymmetric_multichain(str(pdb), keep_chains={'A'}) is False


def test_chain_residue_sequences_dedupes_repeated_atoms_of_one_residue(tmp_path):
    pdb = tmp_path / "in.pdb"
    pdb.write_text(
        _residue_line('A', 1, 'ALA', serial=1) +
        _residue_line('A', 1, 'ALA', serial=2)  # second atom, same residue
    )
    seqs = chain_residue_sequences(str(pdb))
    assert seqs == {'A': ['ALA']}
