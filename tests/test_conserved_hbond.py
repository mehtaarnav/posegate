# posegate/tests/test_conserved_hbond.py
"""Unit tests for multi-residue conserved-contact matching.

find_conserved_hbond accepts a precomputed interaction fingerprint, so
these tests inject a synthetic one rather than building real complexes:
the logic under test is residue matching and any/all mode, not ProLIF's
interaction detection.
"""

import pytest

from posegate.autopsy import find_conserved_hbond, normalize_conserved_residues


class FakeResidue:
    """Stands in for ProLIF's ResidueId (name/number/chain are all it needs)."""

    def __init__(self, name, number, chain):
        self.name = name
        self.number = number
        self.chain = chain


def make_ifp(*residues_with_interactions):
    """Builds a synthetic IFP: each argument is (name, number, chain, kind),
    where kind is 'HBDonor', 'HBAcceptor' or None for a non-hbond contact."""
    ifp = {}
    for i, (name, number, chain, kind) in enumerate(residues_with_interactions):
        meta = {'distance': 2.9, 'DHA_angle': 165.0}
        interactions = {kind: (meta,)} if kind else {'Hydrophobic': ({},)}
        ifp[(FakeResidue('LIG', 1, 'X'), FakeResidue(name, number, chain))] = interactions
    return ifp


def test_normalize_accepts_single_tuple_and_sequence():
    assert normalize_conserved_residues(('ASN', 140, 'A')) == [('ASN', 140, 'A')]
    assert normalize_conserved_residues([('ASN', 140, 'A')]) == [('ASN', 140, 'A')]
    assert normalize_conserved_residues(
        [('GLU', 353, 'A'), ('ARG', 394, 'A')]
    ) == [('GLU', 353, 'A'), ('ARG', 394, 'A')]


def test_normalize_disambiguates_three_residue_sequence():
    """A 3-element sequence of residues must not be mistaken for one
    (name, number, chain) tuple, which is the ambiguous case."""
    three = [('ASP', 25, 'A'), ('ASP', 25, 'B'), ('THR', 199, 'A')]
    assert normalize_conserved_residues(three) == three


def test_single_residue_match():
    ifp = make_ifp(('ASN', 140, 'A', 'HBDonor'))
    hits = find_conserved_hbond(None, None, residues=('ASN', 140, 'A'), ifp=ifp)
    assert len(hits) == 1
    assert hits[0]['residue'] == 'ASN140.A'


def test_chain_is_discriminating():
    """HIV-1 protease's dyad is Asp25 on each chain; matching must not
    collapse the two chains together."""
    ifp = make_ifp(('ASP', 25, 'B', 'HBDonor'))
    assert find_conserved_hbond(None, None, residues=('ASP', 25, 'A'), ifp=ifp) == []
    assert len(find_conserved_hbond(None, None, residues=('ASP', 25, 'B'), ifp=ifp)) == 1


def test_any_mode_satisfied_by_one_of_two():
    ifp = make_ifp(('GLU', 353, 'A', 'HBAcceptor'), ('ARG', 394, 'A', None))
    hits = find_conserved_hbond(
        None, None, residues=[('GLU', 353, 'A'), ('ARG', 394, 'A')], mode='any', ifp=ifp
    )
    assert len(hits) == 1
    assert hits[0]['residue'] == 'GLU353.A'


def test_all_mode_requires_every_residue():
    partial = make_ifp(('GLU', 353, 'A', 'HBAcceptor'), ('ARG', 394, 'A', None))
    assert find_conserved_hbond(
        None, None, residues=[('GLU', 353, 'A'), ('ARG', 394, 'A')], mode='all', ifp=partial
    ) == []

    complete = make_ifp(('GLU', 353, 'A', 'HBAcceptor'), ('ARG', 394, 'A', 'HBDonor'))
    hits = find_conserved_hbond(
        None, None, residues=[('GLU', 353, 'A'), ('ARG', 394, 'A')], mode='all', ifp=complete
    )
    assert {h['residue'] for h in hits} == {'GLU353.A', 'ARG394.A'}


def test_defaults_to_brd4_asn140():
    ifp = make_ifp(('ASN', 140, 'A', 'HBDonor'))
    assert len(find_conserved_hbond(None, None, ifp=ifp)) == 1


def test_rejects_unknown_mode():
    with pytest.raises(ValueError):
        find_conserved_hbond(None, None, residues=('ASN', 140, 'A'), mode='either', ifp={})
