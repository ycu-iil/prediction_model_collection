import pickle
import random

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def plot_hist(pred):
    fig, ax = plt.subplots(figsize=(6,4))
    ax.hist(pred, bins=50)
    ax.set_ylabel("Number of data", fontsize=20, labelpad=20, weight='bold')
    ax.set_xlabel("Predicted value", fontsize=20, labelpad=20, weight='bold')
    ax.xaxis.set_tick_params(direction="out", labelsize=18, width=3, pad=10)
    ax.yaxis.set_tick_params(direction="out", labelsize=18, width=3, pad=10)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.grid()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    return fig


def main():

    # Data preparation
    zinc = Chem.SmilesMolSupplier('../../data/250k_rndm_zinc_drugs_clean.smi')
    zinc_mols = [x for x in zinc]
    random.seed(0)
    mols = random.sample(zinc_mols, 10000)
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
    X = np.array(fps)

    # Load model
    with open('./model/lgb_vegfr2.pkl', mode='rb') as f:
        model = pickle.load(f)

    # Predict
    pred = model.predict(X)

    # Plot histogram
    fig = plot_hist(pred)
    fig.savefig("img/pred_distribution.png", dpi=300, bbox_inches="tight", pad_inches=0.05)


if __name__ == '__main__':
    main()