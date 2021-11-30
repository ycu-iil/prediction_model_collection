import pickle

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def main():

    # Data preparation
    suppl = Chem.SDMolSupplier('./data/metabolic_stability_ChEMBL.sdf')
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in suppl]
    values = []
    for mol in suppl:
        value = mol.GetProp('Metabolic stability (liver micros.- human - 60min)')
        values.append(float(value))
    X = np.array(fps)
    y = np.array(values)

    # Load model
    with open('./model/lgb_hlm.pkl', mode='rb') as f:
        model = pickle.load(f)

    # Predict
    y_pred = model.predict(X)


if __name__ == '__main__':
    main()