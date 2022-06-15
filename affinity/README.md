# Prediction models for inhibitory activities

## Contents

LightGBM regression models that predict inhibitory activity for the following 12 proteins:
- Tyrosine-protein kinase ABL (ABL)
- $\beta$-secretase 1 (BACE1)
- Epidermal growth factor receptor (EGFR)
- Ephrin type-B receptor 4 (EPHB4)
- Receptor protein-tyrosine kinase erbB-2 (ERBB2)
- Receptor protein-tyrosine kinase erbB-4 (ERBB4)
- Fibroblast growth factor receptor 1 (FGFR1)
- Tyrosine-protein kinase HCK (HCK)
- Tyrosine-protein kinase LCK (LCK)
- Platelet-derived growth factor receptor beta (PDGFRbeta)
- Tyrosine-protein kinase SRC (SRC)
- Vascular endothelial growth factor receptor 2 (VEGFR2)

## Datasets

Training data for each model is placed in `./data/`.

### `desc_canvas_aug30.csv`

- Content: Molecules with SMILES and their pIC50 value to BACE1.

- Source: https://github.com/deepchem/deepchem/tree/master/datasets

### `ChEMBL_EGFR.sdf`

- Content: Molecules with SMILES and their pChEMBL value, which is equivalent to the pIC50, to EGFR.

- Source: ChEMBL 28

### `ChEMBL_EGFR_selectivity_data.sdf`

- Content: Molecules with SMILES and their pChEMBL value to the following 11 proteins: ABL, EPHB4, ERBB2, ERBB4, FGFR1, HCK, LCK, PDGFRbeta, SRC, and VEGFR2.

- Source: ChEMBL 28