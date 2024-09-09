import argparse
import os
import pickle
import time

import numpy as np
import optuna
import lightgbm as lgb
from lightgbm import LGBMRegressor
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold, train_test_split


def get_parser():
    parser = argparse.ArgumentParser(
        usage=f'ex.) python {os.path.basename(__file__)}'
    )
    parser.add_argument(
        '--data', type=str, required=True,
        )
    parser.add_argument(
        '--outdir', type=str, required=True,
        )
    parser.add_argument(
        '--n_jobs', type=int, required=False, default=1,
        )
    parser.add_argument(
        '--n_trials', type=int, required=False, default=200,
        )
    parser.add_argument(
        '--debug', action='store_true'
        )
    return parser.parse_args()


def calc_corrcoef(y, y_pred):
    return np.corrcoef(y, y_pred)[0][1]

def calculate_scores(X, y, cv, best_params):
    model_eval = LGBMRegressor(boosting_type='gbdt', objective='regression',
                          random_state=0, n_estimators=100)   # Initialize model
    model_eval.set_params(**best_params)

    r = []
    r2 = []
    mae = []
    mse = []
    rmse = []
    y_train_true_list = []
    y_train_pred_list = []
    y_test_true_list = []
    y_test_pred_list = []

    for (train_val_indices, test_indices) in cv.split(X, y):
        train_indices, val_indices = train_test_split(train_val_indices, test_size=0.25, random_state=0, shuffle=True)
        model_eval.fit(X[train_indices], y[train_indices],
                  eval_metric='mean_squared_error',
                  eval_set=[(X[val_indices], y[val_indices])],
                  callbacks=[lgb.early_stopping(stopping_rounds=10, verbose=True)])
        y_test_pred = model_eval.predict(X[test_indices])
        r.append(calc_corrcoef(y[test_indices], y_test_pred))
        r2.append(r2_score(y[test_indices], y_test_pred))
        mae.append(mean_absolute_error(y[test_indices], y_test_pred))
        mse.append(mean_squared_error(y[test_indices], y_test_pred))
        rmse.append(np.sqrt(mean_squared_error(y[test_indices], y_test_pred)))
        y_train_true_list.append(list(y[train_indices]))
        y_train_pred_list.append(list(model_eval.predict(X[train_indices])))
        y_test_true_list.append(list(y[test_indices]))
        y_test_pred_list.append(list(y_test_pred))

    df_scores = pd.DataFrame({
        'Corr_Coef': [np.mean(r), np.std(r)],
        'R2': [np.mean(r2), np.std(r2)],
        'MAE': [np.mean(mae), np.std(mae)],
        'MSE': [np.mean(mse), np.std(mse)],
        'RMSE': [np.mean(rmse), np.std(rmse)]
        }, index = ['mean', 'std'])

    return df_scores, y_train_true_list, y_train_pred_list, y_test_true_list, y_test_pred_list


def main():

    start_time = time.perf_counter()
    args = get_parser()
    data_path = args.data
    outdir = args.outdir
    n_jobs = args.n_jobs
    n_trials = args.n_trials
    debug = args.debug

    # Data preparation
    df = pd.read_csv(data_path)
    print(f'[INFO] Loaded {data_path} ({len(df)} data points)')
    if debug:
        if len(df) > 100:
            df = df.sample(100, random_state=0)
    mols = [Chem.MolFromSmiles(smi) for smi in df['smiles']]
    values = df['value'].values
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
    descriptors = np.array(fps)
    X = np.array(descriptors)
    y = np.array(values)
    cv = KFold(n_splits=5, shuffle=True, random_state=0)

    # Hyperparameter optimization
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
        scores = []
        for (train_val_indices, test_indices) in cv.split(X, y):   # train:val:test = 7:1:2
            train_indices, val_indices = train_test_split(train_val_indices, test_size=1/8, random_state=0, shuffle=True)
            model.fit(X[train_indices], y[train_indices],
                        eval_metric='mean_squared_error',
                        eval_set=[(X[val_indices], y[val_indices])],
                        callbacks=[lgb.early_stopping(stopping_rounds=10, verbose=True)])
            y_pred = model.predict(X[test_indices])
            score = mean_squared_error(y[test_indices], y_pred)
            scores.append(score)
        return np.mean(scores)

    model = LGBMRegressor(boosting_type='gbdt', objective='regression',
                    random_state=0, n_estimators=100)
    study = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler(seed=0))
    n_trials = 10 if debug else n_trials
    study.optimize(objective, n_trials=n_trials, n_jobs=n_jobs)
    best_params = study.best_trial.params

    # Model fitting (train:val=8:2)
    model_best = LGBMRegressor(boosting_type='gbdt', objective='regression', random_state=0, n_estimators=100)   # Initialize model
    model_best.set_params(**best_params)
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=0, shuffle=True)
    model_best.fit(X_train, y_train,
                eval_metric='mean_squared_error',
                eval_set=[(X_val, y_val)],
                callbacks=[lgb.early_stopping(stopping_rounds=10, verbose=True)])

    # Save model
    os.makedirs(outdir, exist_ok=True)
    model_output_path = os.path.join(outdir, 'lgb.pkl')
    with open(model_output_path, mode='wb') as f:
        pickle.dump(model_best, f)
    print(f'[INFO] model saved at {model_output_path}')

    # Model evaluation with best hyperparameters
    accuracy, y_train_true_list, y_train_pred_list, y_test_true_list, y_test_pred_list = calculate_scores(X, y, cv, best_params)
    accuracy_out_path = f'{outdir}/performance.csv'
    accuracy.to_csv(accuracy_out_path)
    print(f'[INFO] Performance:')
    print(accuracy)
    print(f'[INFO] model performance saved at {accuracy_out_path}')

    list_to_save = [
        ('y_train_true_list.pkl', y_train_true_list),
        ('y_train_pred_list.pkl', y_train_pred_list),
        ('y_test_true_list.pkl', y_test_true_list),
        ('y_test_pred_list.pkl', y_test_pred_list)
        ]
    for filename, l in list_to_save:
        with open(f'{outdir}/{filename}', 'wb') as f:
            pickle.dump(l, f)

    finish_time = time.perf_counter()
    elapsed_time = finish_time - start_time
    elapsed_minute = int(elapsed_time // 60)
    elapsed_second = elapsed_time % 60
    print(f'[INFO] Elapsed time: {elapsed_minute} min {elapsed_second:.1f} sec')
    print('[INFO] Finished!')


if __name__ == '__main__':
    main()