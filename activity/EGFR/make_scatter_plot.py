import os
import pickle
from decimal import Decimal, ROUND_HALF_UP

from lightgbm import LGBMRegressor
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import cross_val_predict, KFold


def make_scatter_plot(y_true, y_pred):
    max_val = max(max(y_true, y_pred, key=max))
    min_val = min(min(y_true, y_pred, key=min))
    corrcoef = np.corrcoef(y_true, y_pred)[0][1]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_pred, y_true, marker='o', s=25, c='dimgray', alpha=0.2)
    ax.set_xlabel('Predicted pChEMBL value', fontsize=20, labelpad=20, weight='normal')
    ax.set_ylabel('Experimental pChEMBL value', fontsize=20, labelpad=20, weight='normal')
    ax.xaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.yaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
    ax.set_xlim(min_val-0.25, max_val+0.25)
    ax.set_ylim(min_val-0.25, max_val+0.25)
    ax.set_aspect('equal')
    ax.grid()
    ax.text(0.8, 0.15,
        f"R = {Decimal(str(corrcoef)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)}",
        transform=ax.transAxes, ha='center', va='center', fontsize=20,
        bbox=(dict(boxstyle='square', fc='white', ec='black', lw='1')))
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    return fig


def main():
    # Data preparation
    suppl = Chem.SDMolSupplier('./data/ChEMBL_EGFR.sdf')
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in suppl]
    values = []
    for mol in suppl:
        value = mol.GetProp('pCV_Epidermal growth factor receptor erbB1')
        values.append(float(value))
    X = np.array(fps)
    y = np.array(values)

    # Load best parameters
    with open('best_parameters.pkl', mode='rb') as f:
        best_params = pickle.load(f)

    # Model creation with best parameters
    seed = 0
    model = LGBMRegressor(boosting_type='gbdt', objective='regression',
                    random_state=seed, n_estimators=100, **best_params)
    fit_params = {
        'verbose': 0,
        'early_stopping_rounds': 10,
        'eval_metric': 'mean_squared_error',
        'eval_set': [(X, y)]
        }

    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    y_pred = cross_val_predict(model, X, y, cv=cv, fit_params=fit_params, method='predict')

    fig = make_scatter_plot(y, y_pred)
    output_path = 'img/'
    os.makedirs(output_path, exist_ok=True)
    image_output_path = os.path.join(output_path, 'scatter_plot.png')
    fig.savefig(image_output_path, dpi=300, bbox_inches='tight', pad_inches=0.05)


if __name__ == '__main__':
    main()