# Epidermal Growth Factor Receptor

## Task

- Regression

- Given a Morgan fingerprint(r=2, 2048 dim), predict the pChEMBL value to EGFR.

## Dataset

- Data size: 6103

<div align="left">
    <img src="img/data_distribution.png" width="400">
</div>

## Model

- LightGBM regressor

- Hyperparameters were optimized in 5-folds cross-validation with Optuna.

- To train the model, run `train.py`.
    - Example usage
        ```bash
        python train.py -o lgb_egfr
        ```

## Accuracy

|Corr Coef|R2|MAE|MSE|RMSE|
|:----:|:----:|:----:|:----:|:----:|
|0.86|0.73|0.49|0.44|0.44|

<div align="left">
      <img src="img/scatter_plot.png" width="400">
</div>