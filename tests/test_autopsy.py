# posegate/tests/test_autopsy.py
import pytest
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from posegate.autopsy import generate_autopsy_report

def test_generate_autopsy_report_integration(tmp_path):
    """Integration test: builds a dummy complex and runs the full autopsy."""
    # Create dummy ligand (benzene)
    lig = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1'))
    AllChem.EmbedMolecule(lig)
    lig_path = tmp_path / "lig.sdf"
    Chem.MolToMolFile(lig, str(lig_path))

    # Create dummy receptor (methane, placed far away to avoid clashes)
    rec = Chem.AddHs(Chem.MolFromSmiles('C'))
    AllChem.EmbedMolecule(rec)
    # Move receptor far away
    conf = rec.GetConformer()
    for i in range(rec.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x + 50, pos.y, pos.z))
        
    rec_path = tmp_path / "rec.pdb"
    Chem.MolToPDBFile(rec, str(rec_path))

    report = generate_autopsy_report(str(lig_path), str(rec_path), vina_score=-5.0)

    assert report['vina_score'] == -5.0
    assert isinstance(report['clashes'], list)
    assert isinstance(report['hbonds'], list)
    assert report['decision'] in ['REJECT', 'PRIORITIZE', 'REVIEW']