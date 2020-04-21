from flask import Flask, request, jsonify
import json
import sys
import requests
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops
from rdkit.Chem import AllChem
import torch
import time
import os

app = Flask(__name__)

os.environ["DGLBACKEND"] = "pytorch"

import dgl
from dgllife.model.model_zoo.gcn_predictor import GCNPredictor
from dgllife.model.model_zoo.acnn import ACNN
from dgllife.model import load_pretrained
from dgllife.utils import smiles_to_bigraph, CanonicalAtomFeaturizer, CanonicalBondFeaturizer, ACNN_graph_construction_and_featurization, get_mol_3d_coordinates

binding_affinity_model = ACNN()
binding_affinity_checkpoint = torch.load('./model/COVID19_binding_affinity.pth')
binding_affinity_model.load_state_dict(binding_affinity_checkpoint)
binding_affinity_model.eval()
print('binding affinity model loaded')
sys.stdout.flush()

toxicity_model = GCNPredictor(in_feats=74,
                             hidden_feats=[64, 64],
                             classifier_hidden_feats=64,
                             n_tasks=12)
toxicity_checkpoint = torch.load('./model/GCN_Tox21_pre_trained.pth')
toxicity_model.load_state_dict(toxicity_checkpoint['model_state_dict'])
toxicity_model.eval()
print('toxicity model loaded')
sys.stdout.flush()


def binding_affinity_pred(ligand_smiles, protein_smiles, model):

    ligand_mol = Chem.MolFromSmiles(ligand_smiles)
    protein_mol = Chem.MolFromSmiles(protein_smiles)

    max_num_ligand_atoms=ligand_mol.GetNumAtoms()
    max_num_protein_atoms=protein_mol.GetNumAtoms()

    AllChem.EmbedMolecule(ligand_mol)
    AllChem.MMFFOptimizeMolecule(ligand_mol)

    AllChem.EmbedMolecule(protein_mol)
    AllChem.MMFFOptimizeMolecule(protein_mol)

    ligand_coordinates = get_mol_3d_coordinates(ligand_mol)
    protein_coordinates = get_mol_3d_coordinates(protein_mol)
    g = ACNN_graph_construction_and_featurization(ligand_mol=ligand_mol, protein_mol=protein_mol, ligand_coordinates=ligand_coordinates , protein_coordinates=protein_coordinates, max_num_ligand_atoms=max_num_ligand_atoms, max_num_protein_atoms=max_num_protein_atoms)

    binding_affinity = float(model(g)[0][0])
    binding_affinity_data = {"binding_affinity" : binding_affinity}
    return binding_affinity_data


def toxicity_pred(smiles, model):
    g = smiles_to_bigraph(smiles ,node_featurizer=CanonicalAtomFeaturizer(), edge_featurizer=CanonicalBondFeaturizer())
    feats = g.ndata['h']
    predictions = model(g, feats)
    predictions = torch.sigmoid(predictions)
    toxicity_data = {"estrogen receptor alpha, LBD (ER, LBD)": None, "estrogen receptor alpha, full (ER, full)": None, "aromatase": None, "aryl hydrocarbon receptor (AhR)": None, "androgen receptor, full (AR, full)": None, "androgen receptor, LBD (AR, LBD)": None, "peroxisome proliferator-activated receptor gamma (PPAR-gamma)": None, "nuclear factor (erythroid-derived 2)-like 2/antioxidant responsive element (Nrf2/ARE)": None, "heat shock factor response element (HSE)": None, "ATAD5": None, "mitochondrial membrane potential (MMP)": None, "p53": None}
    for pred, tox in enumerate(toxicity_data.keys()):
        toxicity_data[tox] = float(predictions[0][pred])
    return toxicity_data


@app.route("/wake", methods=['POST'])
def wake():
    return jsonify({'wake': 'wake'})


@app.route("/binding_affinity", methods=['POST'])
def binding_affinity():
    if request.method == 'POST':
        req_data = request.json
        print(req_data)
        sys.stdout.flush()
        ligand_smiles = req_data['ligand_smiles']
        protein_smiles = req_data['protein_smiles']
        binding_affinity_data = binding_affinity_pred(ligand_smiles, protein_smiles, binding_affinity_model)
        print(binding_affinity_data)
        sys.stdout.flush()
        return jsonify(binding_affinity_data)


@app.route("/toxicity", methods=['POST'])
def toxicity():
    if request.method == 'POST':
        req_data = request.json
        print(req_data)
        sys.stdout.flush()
        smiles = req_data['smiles']
        toxicity_data = toxicity_pred(smiles, toxicity_model)
        print(toxicity_data)
        sys.stdout.flush()
        return jsonify(toxicity_data)


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
