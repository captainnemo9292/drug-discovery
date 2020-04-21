from flask import Flask, request, jsonify
import json
import sys
import requests
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops
import torch
import time
import os

app = Flask(__name__)

os.environ["DGLBACKEND"] = "pytorch"

from dgllife.model.model_zoo.dgmg import DGMG
from dgllife.model.model_zoo.dgmg import MoleculeEnv

atom_types=['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br', 'Mg', 'Na', 'Ca', 'Fe', 'As', 'Al', 'I', 'B', 'V', 'K', 'Tl', 'Yb', 'Sb', 'Sn', 'Ag', 'Pd', 'Co', 'Se', 'Ti', 'Zn', 'H', 'Li', 'Ge', 'Cu', 'Au', 'Ni', 'Cd', 'In', 'Mn', 'Zr', 'Cr', 'Pt', 'Hg', 'Pb']
bond_types=[Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]
node_hidden_size=128
num_prop_rounds=2
dropout=0.2


model = DGMG(atom_types=atom_types,
             bond_types=bond_types,
             node_hidden_size=node_hidden_size,
             num_prop_rounds=num_prop_rounds,
             dropout=dropout)

checkpoint = torch.load('./model/COVID19_molecule_generation.pth')
model.load_state_dict(checkpoint)
print('loaded')
sys.stdout.flush()

def generate_molecule(n_sample, model):
    compound_data = None
    """
    if base_molecule != None:
        print("converting...")
        sys.stdout.flush()

        molecule_env = MoleculeEnv(atom_types=atom_types, bond_types=bond_types)
        base_mol = Chem.MolFromSmiles(base_molecule)
        atom_order = rdmolfiles.CanonicalRankAtoms(base_mol)
        base_mol = rdmolops.RenumberAtoms(base_mol, atom_order)
        actions = molecule_env.get_decision_sequence(base_mol, atom_order)
        print(actions)
        sys.stdout.flush()
    else:
        actions = None
    """
    while compound_data == None:
        generated = model(actions=None, rdkit_mol=True, max_num_steps=400, compute_log_prob=True)
        smiles = generated[1]
        mol = Chem.MolFromSmiles(smiles)
        if mol != None and -generated[0].item()<30:
            print(smiles)
            sys.stdout.flush()
            res = requests.post("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{}/JSON?MaxRecords={}" .format(smiles, n_sample))
            time.sleep(10)
            res = requests.post("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,IUPACName,XLogP,ExactMass,MonoisotopicMass,TPSA,Complexity,Charge,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,IsotopeAtomCount,AtomStereoCount,DefinedAtomStereoCount,UndefinedAtomStereoCount,BondStereoCount,DefinedBondStereoCount,UndefinedBondStereoCount,CovalentUnitCount,Volume3D,XStericQuadrupole3D,YStericQuadrupole3D,ZStericQuadrupole3D,FeatureCount3D,FeatureAcceptorCount3D,FeatureDonorCount3D,FeatureAnionCount3D,FeatureCationCount3D,FeatureRingCount3D,FeatureHydrophobeCount3D,ConformerModelRMSD3D,EffectiveRotorCount3D,ConformerCount3D,Fingerprint2D/JSON" .format(res.json()["Waiting"]["ListKey"]))
            compound_data = res.json()
            print(compound_data)
            sys.stdout.flush()
            if res.status_code != 200:
                compound_data = "failed"
        else:
            pass
    return compound_data


@app.route("/wake", methods=['POST'])
def wake():
    return jsonify({'wake': 'wake'})


@app.route("/generate", methods=['POST'])
def generate():
    if request.method == 'POST':
        req_data = request.json
        print(req_data)
        sys.stdout.flush()
        n_sample = req_data['n_sample']
        compound_data = generate_molecule(n_sample, model)
        print(compound_data)
        sys.stdout.flush()
        return jsonify({'compound_data': compound_data})

if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
