import pickle
from pprint import pprint

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def main():

    # SMILES strings were obtained from https://github.com/molecule-generator-collection/ChemTSv2/blob/master/data/250k_rndm_zinc_drugs_clean.smi
    smiles_list = [
        'CC(C)(C)c1ccc2occ(CC(=O)Nc3ccccc3F)c2c1',
        'C[C@@H]1CC(Nc2cncc(-c3nncn3C)c2)C[C@@H](C)C1',
        'N#Cc1ccc(-c2ccc(O[C@@H](C(=O)N3CCCC3)c3ccccc3)cc2)cc1',
        'CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1',
        'N#CC1=C(SCC(=O)Nc2cccc(Cl)c2)N=C([O-])[C@H](C#N)C12CCCCC2',
        'CC[NH+](CC)[C@](C)(CC)[C@H](O)c1cscc1Br',
        'COc1ccc(C(=O)N(C)[C@@H](C)C/C(N)=N/O)cc1O',
        'O=C(Nc1nc[nH]n1)c1cccnc1Nc1cccc(F)c1',
        'Cc1c(/C=N/c2cc(Br)ccn2)c(O)n2c(nc3ccccc32)c1C#N',
        'C[C@@H]1CN(C(=O)c2cc(Br)cn2C)CC[C@H]1[NH3+]'
        ]
    mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    fp_list = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mol_list]
    X = np.array(fp_list)

    model_path = 'EGFR/lgb.pkl'
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    pred = model.predict(X)

    print('Input SMILES:')
    pprint(smiles_list)
    print('Predicted values:')
    pprint(pred)

if __name__ == '__main__':
    main()