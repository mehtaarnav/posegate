import json
import numpy as np
from rdkit import Chem
from posegate.autopsy import build_ifp
from posegate.receptor_prep import load_receptor_mol

with open('data/five_family/hiv_prepped.json') as f:
    structures = json.load(f)

for s in structures:
    lig = Chem.MolFromMolFile(s['ligand_sdf'], removeHs=False)
    rec = load_receptor_mol(s['receptor_pdb'])
    if lig is None or rec is None:
        print(s['pdb_id'], 'FAILED TO LOAD')
        continue

    ifp = build_ifp(lig, rec)
    for (lres, pres), interactions in ifp.items():
        if pres.number == 25:
            print(f"{s['pdb_id']}: {pres} -> {list(interactions.keys())}")
