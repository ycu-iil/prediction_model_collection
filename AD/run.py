import argparse
import os
import sys
import yaml

import numpy as np
from mordred import Calculator, descriptors
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.preprocessing import StandardScaler

from categorizer import tanimoto_similarity_based, euclidean_distance_based, cosine_similarity_based


def get_parser():
    parser = argparse.ArgumentParser(
        usage=f'ex.) python {os.path.basename(__file__)} -c setting_AD.yaml'
    )
    parser.add_argument(
        '-c', '--config', type=str, required=True,
        help='path to a config file'
        )
    return parser.parse_args()


def standard_scaler(array):
    scaler = StandardScaler()
    scaled_array = scaler.fit_transform(array)
    return scaled_array


def descriptors_calculator(mols):
    calc = Calculator([
        descriptors.Aromatic,
        descriptors.AtomCount,
        descriptors.CarbonTypes.FractionCSP3,
        descriptors.HydrogenBond,
        descriptors.Lipinski,
        descriptors.RingCount.RingCount(aromatic=False, hetero=False),
        descriptors.RingCount.RingCount(aromatic=True),
        descriptors.RingCount.RingCount(hetero=True),
        descriptors.SLogP,
        descriptors.TopoPSA,
        descriptors.RotatableBond,
        descriptors.Weight
        ], ignore_3D=True)
    return calc.pandas(mols).values


def main():
    args = get_parser()
    with open(args.config, "r") as f:
        conf = yaml.load(f, Loader=yaml.SafeLoader)

    # Loading the input data
    input_path = conf['input']['path']
    _, input_ext = os.path.splitext(input_path)
    if input_ext == '.sdf':
        mols = Chem.SDMolSupplier(input_path)
    elif input_ext == '.smi':
        mols = Chem.SmilesMolSupplier(input_path)
    elif input_ext == '.csv':
        df = pd.read_csv(input_path, sep=',')
        mols = [Chem.MolFromSmiles(smi) for smi in df[conf['input']['col_name']]]
    else:
        sys.exit('The file format submitted as input is not supported. Under construction.')
    print(f'[INFO] `{input_path}` has loaded as input.')

    # Loading the reference data
    ref_path = conf['reference']['path']
    _, ref_ext = os.path.splitext(ref_path)
    if ref_ext == '.sdf':
        ref_mols = Chem.SDMolSupplier(ref_path)
    elif ref_ext == '.smi':
        ref_mols = Chem.SmilesMolSupplier(ref_path)
    elif ref_ext == '.csv':
        ref_df = pd.read_csv(ref_path, sep=',')
        ref_mols = [Chem.MolFromSmiles(smi) for smi in ref_df[conf['reference']['col_name']]]
    else:
        sys.exit('The file format submitted as reference is not supported. Under construction.')
    print(f'[INFO] `{ref_path}` has loaded as reference.')

    # Calculation
    ad_params = conf['AD']
    if ad_params['type'] == 'Tanimoto_simirarity':
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in mols]
        ref_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048) for mol in ref_mols]
        result = [tanimoto_similarity_based(fp, ref_fps, conf) for fp in fps]

    elif ad_params['type'] == 'Euclidean_distance-based':
        mrd_raw = descriptors_calculator(mols)
        mrd = standard_scaler(mrd_raw)
        ref_mrd_raw = descriptors_calculator(ref_mols)
        ref_mrd = standard_scaler(ref_mrd_raw)
        result = [euclidean_distance_based(m, ref_mrd, conf) for m in mrd]

    elif ad_params['type'] == 'Cosine_similarity':
        mrd_raw = descriptors_calculator(mols)
        mrd = standard_scaler(mrd_raw)
        ref_mrd_raw = descriptors_calculator(ref_mols)
        ref_mrd = standard_scaler(ref_mrd_raw)
        result = [cosine_similarity_based(m, ref_mrd, conf) for m in mrd]

    else:
        sys.exit('Not supported for AD types other than "Tanimoto_simirarity", "Euclidean_distance-based", and "Cosine_similarity".')

    # Output results to DataFrame
    smiles = [Chem.MolToSmiles(mol) for mol in mols]
    df_result = pd.DataFrame({
        'Smiles': smiles,
        'Result': result
        })

    # Save results
    os.makedirs(conf['output_dir'], exist_ok=True)
    output_path = os.path.join(conf['output_dir'], conf['output_file_name']) + '.csv'
    df_result.to_csv(output_path)
    print(f'[INFO] save results at `{output_path}`')


if __name__ == '__main__':
    main()