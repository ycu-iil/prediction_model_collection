# prediction_model_collection

- This repository is a collection of prediction models for the inhibitory activity and ADMET properties of molecules.

- For inhibitory activities, 12 prediction models for the following proteins, including 11 tyrosine kinases, are stored.

  - Abelson tyrosine-protein kinase (ABL)
  - Epidermal growth factor receptor (EGFR)
  - Ephrin type-B receptor 4 (EPHB4)
  - Receptor protein-tyrosine kinase erbB-2 (ERBB2)
  - Receptor protein-tyrosine kinase erbB-4 (ERBB4)
  - Fibroblast growth factor receptor 1 (FGFR1)
  - Tyrosine-protein kinase HCK (HCK)
  - Tyrosine-protein kinase Lck (LCK)
  - Platelet-derived growth factor receptor beta (PDGFRbeta)
  - Proto-oncogene tyrosine-protein kinase Src (SRC)
  - Vascular endothelial growth factor receptor 2 (VEGFR2)
  - $\beta$-secretase 1 (BACE1)

- For ADMET properties, 4 prediction models for the following properties, are stored.

  - Caco-2 permeability (Permeability)
  - Metabolic stability in human liver microsomes (Metabolic stability)
  - Solubility
  - Oral rat acute toxicity (Toxicity)


## Dependency

### Environment (confirmed)

- CentOS: 7.8
- Ubuntu: 20.04.2 LTS
- macOS Monterey: 12.4

## Package

1. python: 3.7
2. RDKit: 2021.03.5
3. PyTDC: 0.3.2
4. LightGBM: 3.2.1
5. Optuna: 2.10.0
6. pandas

## How to setup (example)

```bash
conda create -n pred_models -c conda-forge python=3.7
# switch a python virtual environment to `pred_models`
conda install -c conda-forge rdkit=2021.03.5 pytdc=0.3.2 lightgbm=3.2.1 optuna=2.10.0 pandas
```

## License

This package is distributed under the MIT License.

## Contact

- Tatsuya Yoshizawa (w225435a@yokohama-cu.ac.jp)
- Shoichi Ishida (ishida.sho.nm@yokohama-cu.ac.jp)