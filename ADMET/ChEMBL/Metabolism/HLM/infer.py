import os
import pickle

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def make_scatter_plot(y_true, y_pred):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_pred, y_true, marker='o', s=25, c='dimgray', alpha=0.2)
    ax.set_xlabel('Predicted value', fontsize=20, labelpad=20, weight='bold')
    ax.set_ylabel('Experimental value', fontsize=20, labelpad=20, weight='bold')
    ax.xaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.yaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(25))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(25))
    ax.set_aspect('equal')
    ax.grid()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    return fig


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

    fig = make_scatter_plot(y, y_pred)
    output_path = 'img/'
    os.makedirs(output_path, exist_ok=True)
    image_output_path = os.path.join(output_path, 'scatter_plot.png')
    fig.savefig(image_output_path, dpi=300, bbox_inches='tight', pad_inches=0.05)


if __name__ == '__main__':
    main()