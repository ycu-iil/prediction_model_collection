import pickle

ABL_MODEL_PATH = '../affinity/ABL/model/lgb_abl.pkl'
EGFR_MODEL_PATH = '../affinity/EGFR/model/lgb_egfr.pkl'
EPHB4_MODEL_PATH = '../affinity/EPHB4/model/lgb_ephb4.pkl'
ERBB2_MODEL_PATH = '../affinity/ERBB2/model/lgb_erbb2.pkl'
ERBB4_MODEL_PATH = '../affinity/ERBB4/model/lgb_erbb4.pkl'
FGFR1_MODEL_PATH = '../affinity/FGFR1/model/lgb_fgfr1.pkl'
HCK_MODEL_PATH = '../affinity/HCK/model/lgb_hck.pkl'
LCK_MODEL_PATH = '../affinity/LCK/model/lgb_lck.pkl'
PDGFRBETA_MODEL_PATH = '../affinity/PDGFRbeta/model/lgb_pdgfrbeta.pkl'
SRC_MODEL_PATH = '../affinity/SRC/model/lgb_src.pkl'
VEGFR2_MODEL_PATH = '../affinity/VEGFR2/model/lgb_vegfr2.pkl'

PERMEABILITY_MODEL_PATH = '../ADMET/ChEMBL/Absorption/Caco-2/model/lgb_caco-2.pkl'
METABOLISM_MODEL_PATH = '../ADMET/ChEMBL/Metabolism/HLM/model/lgb_hlm.pkl'
SOLUBILITY_MODEL_PATH = '../ADMET/TDC/Absorption/Solubility/model/lgb_sol.pkl'
TOXICITY_MODEL_PATH = '../ADMET/TDC/Toxicity/Acute_Toxicity_LD50/model/lgb_tox.pkl'


with open(ABL_MODEL_PATH, mode='rb') as abl, \
     open(EGFR_MODEL_PATH, mode='rb') as egfr, \
     open(EPHB4_MODEL_PATH, mode='rb') as ephb4, \
     open(ERBB2_MODEL_PATH, mode='rb') as erbb2, \
     open(ERBB4_MODEL_PATH, mode='rb') as erbb4, \
     open(FGFR1_MODEL_PATH, mode='rb') as fgfr1, \
     open(HCK_MODEL_PATH, mode='rb') as hck, \
     open(LCK_MODEL_PATH, mode='rb') as lck, \
     open(PDGFRBETA_MODEL_PATH, mode='rb') as pdgfrbeta, \
     open(SRC_MODEL_PATH, mode='rb') as src, \
     open(VEGFR2_MODEL_PATH, mode='rb') as vegfr2, \
     open(PERMEABILITY_MODEL_PATH, mode='rb') as perm, \
     open(SOLUBILITY_MODEL_PATH, mode='rb') as sol, \
     open(METABOLISM_MODEL_PATH, mode='rb') as meta, \
     open(TOXICITY_MODEL_PATH, mode='rb') as tox :

    lgb_abl = pickle.load(abl)
    print(f"[INFO] loaded model from {ABL_MODEL_PATH}")
    lgb_egfr = pickle.load(egfr)
    print(f"[INFO] loaded model from {EGFR_MODEL_PATH}")
    lgb_ephb4 = pickle.load(ephb4)
    print(f"[INFO] loaded model from {EPHB4_MODEL_PATH}")
    lgb_erbb2 = pickle.load(erbb2)
    print(f"[INFO] loaded model from {ERBB2_MODEL_PATH}")
    lgb_erbb4 = pickle.load(erbb4)
    print(f"[INFO] loaded model from {ERBB4_MODEL_PATH}")
    lgb_fgfr1 = pickle.load(fgfr1)
    print(f"[INFO] loaded model from {FGFR1_MODEL_PATH}")
    lgb_hck = pickle.load(hck)
    print(f"[INFO] loaded model from {HCK_MODEL_PATH}")
    lgb_lck = pickle.load(lck)
    print(f"[INFO] loaded model from {LCK_MODEL_PATH}")
    lgb_pdgfrbeta = pickle.load(pdgfrbeta)
    print(f"[INFO] loaded model from {PDGFRBETA_MODEL_PATH}")
    lgb_src = pickle.load(src)
    print(f"[INFO] loaded model from {SRC_MODEL_PATH}")
    lgb_vegfr2 = pickle.load(vegfr2)
    print(f"[INFO] loaded model from {VEGFR2_MODEL_PATH}")

    lgb_perm = pickle.load(perm)
    print(f"[INFO] loaded model from {PERMEABILITY_MODEL_PATH}")
    lgb_sol = pickle.load(sol)
    print(f"[INFO] loaded model from {SOLUBILITY_MODEL_PATH}")
    lgb_meta = pickle.load(meta)
    print(f"[INFO] loaded model from {METABOLISM_MODEL_PATH}")
    lgb_tox = pickle.load(tox)
    print(f"[INFO] loaded model from {TOXICITY_MODEL_PATH}")


def main():
    models = {
        'ABL': lgb_abl,
        'EGFR': lgb_egfr,
        'EPHB4': lgb_ephb4,
        'ERBB2': lgb_erbb2,
        'ERBB4': lgb_erbb4,
        'FGFR1': lgb_fgfr1,
        'HCK': lgb_hck,
        'LCK': lgb_lck,
        'PDGFRbeta': lgb_pdgfrbeta,
        'SRC': lgb_src,
        'VEGFR2': lgb_vegfr2,
        'Sol': lgb_sol,
        'Perm': lgb_perm,
        'Meta': lgb_meta,
        'Tox': lgb_tox
        }
    with open("lgb_models.pickle", "wb") as f:
        pickle.dump(models, f)


if __name__ == "__main__":
    main()