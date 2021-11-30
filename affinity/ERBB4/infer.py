import pickle

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def main():

    # Data preparation
    suppl = Chem.SDMolSupplier("./data/ChEMBL_EGFR_selectivity_data.sdf")
    values = []
    mols = []
    for mol in suppl:
        if mol is None:
            continue
        try:
            value = mol.GetProp('pCV_Receptor protein-tyrosine kinase erbB-4')
        except:
            continue
        mols.append(mol)
        values.append(float(value))
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
    X = np.array(fps)
    y = np.array(values)

    # Load model
    with open('./model/lgb_erbb4.pkl', mode='rb') as f:
        model = pickle.load(f)

    # Predict
    y_pred = model.predict(X)


if __name__ == '__main__':
    main()