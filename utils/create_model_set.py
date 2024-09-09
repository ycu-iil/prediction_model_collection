import glob
import pickle


def main():
    model_path = glob.glob('../*/lgb.pkl')
    target_names = [path.split('/')[1] for path in model_path]
    rename_dict = {'Metabolic_stability': 'Stab', 'Solubility': 'Sol', 'Permeability': 'Perm'}
    skip_list = ['BACE1', 'ERBB4', 'HCK']

    model_set = {}
    for t, m in zip(target_names, model_path):
        if t in skip_list:
            continue
        if t in rename_dict:
            t = rename_dict[t]
        with open(m, 'rb') as f:
            model = pickle.load(f)
        model_set[t] = model

    with open('lgb_models.pkl', 'wb') as f:
        pickle.dump(model_set, f)


if __name__ == '__main__':
    main()