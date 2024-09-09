# Prediction Model Collection

This repository conteins a collection of prediction models for some properties of small molecules.


## Packages

- Python: 3.11
- LightGBM: 4.5.0
- RDKit: 2023.9.2
- Optuna: 3.6.1
- NumPy: 1.26.4
- Matplotlib: 3.9.0
- pandas: 2.2.2
- scikit-learn: 1.5.1


## How to run

- `train.py` is a script to build a LightGBM regressor.
- Hyperparameters are optimized in 5-folds cross-validation with Optuna.
- To train the model, run `train.py` like the following.
    ```bash
    python run.py --data ./data/EGFR.csv --outdir ./EGFR
    ```
  - `--data`: Path to data to be used to train the model.
  - `--outdir`: Path to where a trained model and a predictive performance report are stored.
- If you want to predict properties of your own molecule, please refer to `infer.py` for how to use the models.


## License

This package is distributed under the MIT License.


## Contact

- Tatsuya Yoshizawa (tatsuya.yoshizawa@riken.jp)
- Kei Terayama (terayama@yokohama-cu.ac.jp)