# posegate/tests/test_audit_fixes.py
"""Regression tests, one per fix, from a code-quality audit response
(see AUDIT_RESPONSE.md). Each test reproduces the specific failure mode
described there and would fail against the pre-fix code.
"""

import math

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from posegate.docking import require_obabel
from posegate.receptor_prep import pdb_atom_name, find_discarded_heterogens
from posegate.autopsy import find_metal_coordination
from posegate.conserved_contacts import wilson_interval


# --- P1: obabel presence check -------------------------------------------

def test_require_obabel_raises_with_actionable_message(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(RuntimeError, match="obabel"):
        require_obabel()


def test_require_obabel_passes_when_present(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/obabel")
    require_obabel()  # must not raise


# --- P2: atom-name truncation ---------------------------------------------

def test_pdb_atom_name_does_not_truncate_four_char_names():
    """The old formatting, f' {name:<3}'[:4], padded 4-char names with a
    leading space and then cut the result back to 4 characters, dropping
    the last character. HD11 and HD12 both became 'HD1'."""
    assert pdb_atom_name("HD11") == "HD11"
    assert pdb_atom_name("HD12") == "HD12"
    assert pdb_atom_name("HD11") != pdb_atom_name("HD12")


def test_pdb_atom_name_still_pads_short_names():
    assert pdb_atom_name("CA") == " CA "
    assert pdb_atom_name("OG1") == " OG1"


# --- P3: metal coordination (reporting-only) ------------------------------

def _atom_with_position(symbol, xyz, resname=None, resnum=1, chain='A'):
    mol = Chem.RWMol()
    idx = mol.AddAtom(Chem.Atom(symbol))
    if resname:
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f' {symbol:<3}')
        info.SetResidueName(resname)
        info.SetResidueNumber(resnum)
        info.SetChainId(chain)
        mol.GetAtomWithIdx(idx).SetMonomerInfo(info)
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, xyz)
    mol.AddConformer(conf)
    return mol.GetMol()


def test_metal_coordination_detected_within_cutoff():
    receptor = _atom_with_position('Zn', (0.0, 0.0, 0.0), resname='ZN', resnum=263, chain='A')
    ligand = _atom_with_position('N', (2.0, 0.0, 0.0))  # 2.0 A: real coordination distance
    hits = find_metal_coordination(ligand, receptor)
    assert len(hits) == 1
    assert hits[0]['metal'] == 'ZN263.A'
    assert hits[0]['distance_A'] == pytest.approx(2.0, abs=0.01)


def test_metal_coordination_absent_beyond_cutoff():
    receptor = _atom_with_position('Zn', (0.0, 0.0, 0.0), resname='ZN', resnum=263, chain='A')
    ligand = _atom_with_position('N', (4.0, 0.0, 0.0))  # 4.0 A: not coordinating
    assert find_metal_coordination(ligand, receptor) == []


def test_metal_coordination_ignores_non_coordinating_elements():
    """A carbon at coordination distance should not be reported: only
    N/O/S commonly donate a lone pair to a coordinated metal."""
    receptor = _atom_with_position('Zn', (0.0, 0.0, 0.0), resname='ZN', resnum=263, chain='A')
    ligand = _atom_with_position('C', (2.0, 0.0, 0.0))
    assert find_metal_coordination(ligand, receptor) == []


def test_metal_coordination_no_metal_in_receptor():
    receptor = _atom_with_position('C', (0.0, 0.0, 0.0))
    ligand = _atom_with_position('N', (2.0, 0.0, 0.0))
    assert find_metal_coordination(ligand, receptor) == []


# --- P4: silent cofactor deletion ------------------------------------------

def test_find_discarded_heterogens_flags_cofactor(tmp_path):
    """A heme iron's HEM residue must be flagged: PDBFixer's
    removeHeterogens would delete it and receptor_prep does not put
    non-metal cofactors back, unlike the retained metal ions."""
    pdb = tmp_path / "hem.pdb"
    pdb.write_text(
        "HETATM    1  FE  HEM A 200      10.000  10.000  10.000  1.00  0.00          FE\n"
        "HETATM    2  O   HOH A 300      20.000  20.000  20.000  1.00  0.00           O\n"
        "HETATM    3 ZN   ZN  A 262       5.000   5.000   5.000  1.00  0.00          ZN\n"
    )
    discarded = find_discarded_heterogens(str(pdb))
    assert discarded == {"HEM"}  # water and the retained metal are excluded


def test_find_discarded_heterogens_empty_when_nothing_to_warn_about(tmp_path):
    pdb = tmp_path / "clean.pdb"
    pdb.write_text(
        "HETATM    1  O   HOH A 300      20.000  20.000  20.000  1.00  0.00           O\n"
        "HETATM    2 ZN   ZN  A 262       5.000   5.000   5.000  1.00  0.00          ZN\n"
    )
    assert find_discarded_heterogens(str(pdb)) == set()


# --- P5: frequency has no notion of statistical weight ---------------------

def test_wilson_interval_distinguishes_equal_frequencies_at_different_n():
    """3/6 and 11/22 are both a raw frequency of 0.5, but 3/6 is far less
    certain. The interval widths must reflect that, or the miner's
    frequency threshold treats them as equally conclusive."""
    lo6, hi6 = wilson_interval(3, 6)
    lo22, hi22 = wilson_interval(11, 22)
    assert (hi6 - lo6) > (hi22 - lo22)
    # both intervals still center near 0.5
    assert abs((lo6 + hi6) / 2 - 0.5) < 0.05
    assert abs((lo22 + hi22) / 2 - 0.5) < 0.05


def test_wilson_interval_bounds_are_valid_probabilities():
    for count, n in [(0, 6), (6, 6), (3, 6), (1, 100)]:
        lo, hi = wilson_interval(count, n)
        assert 0.0 <= lo <= hi <= 1.0


# --- P6: scaler must be refit per fold/resample, not fit once globally ----

def test_compare_feature_weights_scaler_is_inside_the_pipeline():
    """Regression guard against re-introducing the global-fit leak: the
    scaler must live inside the cross-validated model, not be applied to
    the whole dataset before cross_val_predict/bootstrap resampling see it."""
    import scripts.compare_feature_weights as cfw
    model = cfw.make_model(seed=0)
    assert "standardscaler" in model.named_steps
    assert "logisticregressioncv" in model.named_steps


def test_fit_once_accepts_raw_unscaled_features():
    """fit_once must take raw features and scale internally, not expect a
    caller to have pre-scaled them (that pre-scaling is exactly the leak
    this fix removes)."""
    import numpy as np
    import scripts.compare_feature_weights as cfw

    rng = np.random.default_rng(0)
    n = 40
    X_raw = np.column_stack([
        rng.normal(-7, 1, n),      # vina_score-like: large negative scale
        rng.integers(0, 4, n),     # hbond_count-like: small integer scale
        rng.integers(0, 2, n),     # conserved_hbond-like: binary
        rng.integers(0, 2, n),
        rng.integers(0, 2, n),
    ]).astype(float)
    y = np.array(([0, 1] * (n // 2)))
    coefs = cfw.fit_once(X_raw, y, seed=0)
    assert len(coefs) == 5
    assert np.all(np.isfinite(coefs))


# --- P7: default score weights match the current, documented BRD4 fit -----

def test_compute_posegate_score_signs_match_current_brd4_fit():
    """The hardcoded posegate_score weights are provisional BRD4 defaults
    that must be updated whenever recalibrate_weights.py is rerun on a new
    BRD4 benchmark. This pins their *direction* (not magnitude, which
    would make the test brittle to a legitimate refit): the specific
    conserved contact rewards (lowers) the score, while generic H-bond
    count, clash count and aromatic count all penalize (raise) it in the
    current fit."""
    from posegate.autopsy import compute_posegate_score

    baseline = compute_posegate_score(-5.0, n_hbonds=0, n_clashes=0, n_aromatic=0, conserved_hit=False)
    assert baseline == pytest.approx(-5.0)

    with_conserved = compute_posegate_score(-5.0, n_hbonds=0, n_clashes=0, n_aromatic=0, conserved_hit=True)
    assert with_conserved < baseline  # a reward: lowers the score

    with_hbond = compute_posegate_score(-5.0, n_hbonds=1, n_clashes=0, n_aromatic=0, conserved_hit=False)
    with_clash = compute_posegate_score(-5.0, n_hbonds=0, n_clashes=1, n_aromatic=0, conserved_hit=False)
    with_aromatic = compute_posegate_score(-5.0, n_hbonds=0, n_clashes=0, n_aromatic=1, conserved_hit=False)
    assert with_hbond > baseline
    assert with_clash > baseline
    assert with_aromatic > baseline


def test_compute_posegate_score_uses_current_weight_magnitudes():
    """Pins the exact constants too, so a silent edit to
    POSEGATE_SCORE_WEIGHTS (as opposed to a deliberate, documented refit)
    is caught."""
    from posegate.autopsy import POSEGATE_SCORE_WEIGHTS
    assert POSEGATE_SCORE_WEIGHTS == {
        'hbond_count': 1.047,
        'clash_count': 5.294,
        'aromatic_count': 6.584,
        'conserved_hbond': -2.063,
    }
