# Prediction models for ADMET properties

## Contents

LightGBM regression models that predict the following four ADMET properties:

- Caco-2 permeability
- Metabolic stability in human liver microsomes
- Solubility
- Oral rat acute toxicity

## Directory Structure

The contents are divided into two directories according to the two sources of training data, ChEMBL 28 and Therapeutics Data Commons (TDC).

A simplified directory structure is shown below.

```
├── ChEMBL
│   ├── Permeability
│   └── Metabolic_stability
├── TDC
│   ├── Solubility
│   └── Toxicity
└── README.md
```