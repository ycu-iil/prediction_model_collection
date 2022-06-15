import argparse
import os
import pickle

from lightgbm import LGBMRegressor
import numpy as np
import optuna
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_validate, cross_val_score, KFold


def get_parser():
    parser = argparse.ArgumentParser(
        usage=f'ex.) python {os.path.basename(__file__)} -o lgb_fgfr1'
    )
    parser.add_argument(
        '-o', '--output', type=str, required=False, default='lgb_fgfr1',
        help='file name of the output model'
        )
    return parser.parse_args()


def calculate_scores(model, X, y, cv):
    def calc_corrcoef(y, y_pred):
        return np.corrcoef(y, y_pred)[0][1]
    scoring = {
        'Corr_Coef': make_scorer(calc_corrcoef),
        'r2':'r2',
        'MAE': 'neg_mean_absolute_error',
        'MSE': 'neg_mean_squared_error',
        'RMSE': 'neg_root_mean_squared_error'
        }
    scores = cross_validate(model, X, y,
                            cv=cv, scoring=scoring, n_jobs=-1)
    df_scores_mean = pd.DataFrame({'Corr_Coef': scores['test_Corr_Coef'].mean(),
                          'R2': scores['test_r2'].mean(),
                          'MAE': -1 * scores['test_MAE'].mean(),
                          'MSE': -1 * scores['test_MSE'].mean(),
                          'RMSE': -1 * scores['test_RMSE'].mean()},
                           index = ['scores'])
    return df_scores_mean


def main():
    def objective(trial):
        params = {
            'lambda_l1': trial.suggest_float('lambda_l1', 1e-8, 10.0, log=True),
            'lambda_l2': trial.suggest_float('lambda_l2', 1e-8, 10.0, log=True),
            'num_leaves': trial.suggest_int('num_leaves', 2, 256),
            'feature_fraction': trial.suggest_float('feature_fraction', 0.4, 1.0),
            'bagging_fraction': trial.suggest_float('bagging_fraction', 0.4, 1.0),
            'bagging_freq': trial.suggest_int('bagging_freq', 1, 7),
            'min_child_samples': trial.suggest_int('min_child_samples', 5, 100)
        }
        model.set_params(**params)
        scores = cross_val_score(model, X, y, cv=cv,
                                scoring='neg_mean_squared_error', fit_params=fit_params, n_jobs=-1)
        return scores.mean()

    args = get_parser()

    # Data preparation
    suppl = Chem.SDMolSupplier("./../data/ChEMBL_EGFR_selectivity_data.sdf")
    mols = []
    values = []
    for mol in suppl:
        if mol is None:
            continue
        try:
            value = mol.GetProp('pCV_Fibroblast growth factor receptor 1')
        except:
            continue
        mols.append(mol)
        values.append(float(value))
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
    X = np.array(fps)
    y = np.array(values)

    # Model creation
    seed = 0
    model = LGBMRegressor(boosting_type='gbdt', objective='regression',
                    random_state=seed, n_estimators=100)
    fit_params = {
        'verbose': 0,
        'early_stopping_rounds': 10,
        'eval_metric': 'mean_squared_error',
        'eval_set': [(X, y)]
        }
    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    study = optuna.create_study(direction='maximize',sampler=optuna.samplers.TPESampler(seed=seed))
    study.optimize(objective, n_trials=200)
    best_params = study.best_trial.params
    model.set_params(**best_params)
    model.fit(X, y)

    # Save best parameters
    with open('best_parameters.pkl', mode='wb') as f:
        pickle.dump(best_params, f)

    # Model evaluation
    accuracy = calculate_scores(model, X, y, cv)
    accuracy.to_csv('accuracy.csv')

    # Save model
    model_path = 'model/'
    os.makedirs(model_path, exist_ok=True)
    model_output_path = os.path.join(model_path, f'{args.output}.pkl')
    with open(model_output_path, mode='wb') as f:
        pickle.dump(model, f)


if __name__ == '__main__':
    main()