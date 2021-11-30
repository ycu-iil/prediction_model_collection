import pickle

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def main():

    # Data preparation
    df = pd.read_csv('./data/desc_canvas_aug30.csv', sep=',')
    mols = [Chem.MolFromSmiles(smi) for smi in df['mol']]
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
    X = np.array(fps)
    y = np.array(df['pIC50'])

    # Load model
    with open('./model/lgb_bace1.pkl', mode='rb') as f:
        model = pickle.load(f)

    # Predict
    y_pred = model.predict(X)


if __name__ == '__main__':
    main()