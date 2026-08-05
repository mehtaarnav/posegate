# posegate/tests/test_conserved_contacts.py
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from posegate.conserved_contacts import mine_conserved_contacts, leave_one_out_validate

def _make_structure(tmp_path, name, shift):
    """A tiny synthetic 'structure': a methanol ligand H-bonding to a water
    receptor placed at a fixed offset, so the same (residue, interaction)
    pair should be mined as conserved across repeated calls with the same
    shift, self-contained (no downloaded PDB data required)."""
    methanol = Chem.AddHs(Chem.MolFromSmiles('CO'))
    AllChem.EmbedMolecule(methanol, randomSeed=1)
    donor_o_idx = 1
    donor_h_idx = [n.GetIdx() for n in methanol.GetAtomWithIdx(donor_o_idx).GetNeighbors() if n.GetSymbol() == 'H'][0]

    water = Chem.AddHs(Chem.MolFromSmiles('O'))
    AllChem.EmbedMolecule(water, randomSeed=1)

    conf = methanol.GetConformer()
    import numpy as np
    pos_o = np.array(conf.GetAtomPosition(donor_o_idx))
    pos_h = np.array(conf.GetAtomPosition(donor_h_idx))
    direction = (pos_h - pos_o) / np.linalg.norm(pos_h - pos_o)
    target_pos_a = pos_h + direction * 2.0

    water_conf = water.GetConformer()
    water_o_pos = np.array(water_conf.GetAtomPosition(0))
    move = target_pos_a - water_o_pos + np.array(shift)
    for atom in water.GetAtoms():
        p = water_conf.GetAtomPosition(atom.GetIdx())
        water_conf.SetAtomPosition(atom.GetIdx(), (p.x + move[0], p.y + move[1], p.z + move[2]))

    lig_path = tmp_path / f"{name}_lig.sdf"
    rec_path = tmp_path / f"{name}_rec.pdb"
    Chem.MolToMolFile(methanol, str(lig_path))
    Chem.MolToPDBFile(water, str(rec_path))
    return {'pdb_id': name, 'ligand_sdf': str(lig_path), 'receptor_pdb': str(rec_path)}

def test_mine_conserved_contacts_finds_repeated_interaction(tmp_path):
    structures = [_make_structure(tmp_path, f"s{i}", shift=(0, 0, 0)) for i in range(3)]
    results = mine_conserved_contacts(structures)

    assert len(results) > 0
    top = results[0]
    assert top['frequency'] == 1.0
    assert top['n_structures'] == 3

def test_mine_conserved_contacts_reports_partial_frequency(tmp_path):
    # Two structures with the same water position (should H-bond), one with
    # the water pushed far away (should not) -> that interaction should
    # show up at 2/3 frequency, not 3/3.
    structures = [
        _make_structure(tmp_path, "near1", shift=(0, 0, 0)),
        _make_structure(tmp_path, "near2", shift=(0, 0, 0)),
        _make_structure(tmp_path, "far", shift=(20, 20, 20)),
    ]
    results = mine_conserved_contacts(structures)

    hbond_hits = [r for r in results if 'HB' in r['interaction']]
    assert len(hbond_hits) > 0
    assert hbond_hits[0]['n_structures'] == 2
    assert hbond_hits[0]['frequency'] == pytest.approx(2 / 3, abs=0.01)


def test_leave_one_out_validate_predicts_a_contact_present_in_every_fold(tmp_path):
    """All five structures share the same H-bond geometry, so whichever
    one is held out, the miner trained on the other four must still
    predict a residue the held-out structure actually contacts: top-1
    accuracy should be perfect. This is the ordinary case the feature
    exists for -- a genuinely conserved contact -- checked with no
    ground truth beyond the ensemble itself, exactly as a real user's
    unstudied target would be."""
    structures = [_make_structure(tmp_path, f"s{i}", shift=(0, 0, 0)) for i in range(5)]
    result = leave_one_out_validate(structures, top_k=(1,))

    assert result['n_ensemble'] == 5
    assert result['n_usable'] == 5
    assert result['accuracy'][1]['n_folds'] == 5
    assert result['accuracy'][1]['accuracy'] == 1.0
    for fold in result['folds']:
        assert fold['top1_hit'] is True


def test_leave_one_out_validate_reports_misses_for_a_one_off_contact(tmp_path):
    """One structure's water is pushed far away, so it shares no
    interaction with the other four at all: whichever fold holds out one
    of the four near structures, the far structure is never mined as a
    prediction (it was never conserved across the other three), but more
    importantly, holding out the far structure itself must miss, since
    nothing trained on the four near structures predicts residues at the
    far structure's actual (near-empty) contact set."""
    structures = [
        _make_structure(tmp_path, "near1", shift=(0, 0, 0)),
        _make_structure(tmp_path, "near2", shift=(0, 0, 0)),
        _make_structure(tmp_path, "near3", shift=(0, 0, 0)),
        _make_structure(tmp_path, "near4", shift=(0, 0, 0)),
        _make_structure(tmp_path, "far", shift=(20, 20, 20)),
    ]
    result = leave_one_out_validate(structures, top_k=(1,))

    by_id = {f['pdb_id']: f for f in result['folds']}
    assert by_id['far']['held_out_residues'] == []
    assert by_id['far']['top1_hit'] is False
    # accuracy is 4/5: every 'near' fold hits, the 'far' fold misses.
    assert result['accuracy'][1]['accuracy'] == pytest.approx(4 / 5)


def test_leave_one_out_validate_handles_a_structure_that_fails_to_load(tmp_path):
    structures = [_make_structure(tmp_path, f"s{i}", shift=(0, 0, 0)) for i in range(3)]
    structures.append({'pdb_id': 'broken', 'ligand_sdf': str(tmp_path / 'nope.sdf'),
                       'receptor_pdb': str(tmp_path / 'nope.pdb')})

    result = leave_one_out_validate(structures, top_k=(1,))

    assert result['n_ensemble'] == 4
    assert result['n_usable'] == 3
    broken_fold = next(f for f in result['folds'] if f['pdb_id'] == 'broken')
    assert broken_fold['skipped'] is True
    assert result['accuracy'][1]['n_folds'] == 3
