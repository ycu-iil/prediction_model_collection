# Vascular endothelial growth factor receptor 2

## Task

- Regression

- Given a Morgan fingerprint(r=2, 2048 dim), predict the pChEMBL value to VEGFR2.

## Dataset

- Data size: 8001

<div align="left">
    <img src="img/data_distribution.png" width="400">
</div>

## Model

- LightGBM regressor

- Hyperparameters were optimized in 5-folds cross-validation with Optuna.

- To train the model, run `train.py`.
    - Example usage
        ```bash
        python train.py -o lgb_vegfr2
        ```

## Accuracy

|Corr Coef|R2|MAE|MSE|RMSE|
|:----:|:----:|:----:|:----:|:----:|
|0.81|0.66|0.46|0.39|0.63|

<div align="left">
      <img src="img/scatter_plot.png" width="400">
</div>
