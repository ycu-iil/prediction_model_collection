import pickle

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from tdc.single_pred import Tox


def smiles_to_morganfp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)


def main():

    # Data preparation
    data = Tox(name = 'LD50_Zhu')
    df = data.get_data(format='df')
    df['MorganFP'] = df['Drug'].apply(lambda x: smiles_to_morganfp(x))
    X = np.array(df['MorganFP'].tolist())
    y = np.array(df['Y'])

    # Load model
    with open('./model/lgb_tox.pkl', mode='rb') as f:
        model = pickle.load(f)

    # Predict
    y_pred = model.predict(X)


if __name__ == '__main__':
    main()