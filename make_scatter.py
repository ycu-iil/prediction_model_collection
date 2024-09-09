import argparse
from decimal import Decimal, ROUND_HALF_UP
import os
import pickle
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


def get_parser():
    parser = argparse.ArgumentParser(
        usage=f'ex.) python {os.path.basename(__file__)}'
    )
    parser.add_argument(
        '--result_dir', type=str, required=True,
        )
    parser.add_argument(
        '--format', type=str, default='png',
        )
    return parser.parse_args()


def make_scatter_plot(y_true, y_pred, corrcoef_ave, corrcoef_std):
    max_val = max(max(y_true, y_pred, key=max))
    min_val = min(min(y_true, y_pred, key=min))
    # max_val = 200
    axis_interval = math.ceil((max_val-min_val)/10)
    axis_space = axis_interval / 4
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(y_pred, y_true, marker='o', s=20, c='dimgray', alpha=0.2)
    ax.set_xlabel('Predicted value', fontsize=20, labelpad=15, weight='normal')
    ax.set_ylabel('Experimental value', fontsize=20, labelpad=15, weight='normal')
    ax.xaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.yaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(axis_interval))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(axis_interval))
    ax.set_xlim(min_val-axis_space, max_val+axis_space)
    ax.set_ylim(min_val-axis_space, max_val+axis_space)
    ax.set_aspect('equal')
    ax.grid(True, which='both', linestyle='-', linewidth=2)
    ax.text(0.8, 0.15,
        f"R = {Decimal(str(corrcoef_ave)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)}"
        f" ± {Decimal(str(corrcoef_std)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)}",
        transform=ax.transAxes, ha='center', va='center', fontsize=15, alpha=1,
        bbox=(dict(boxstyle='square', fc='white', ec='gray', lw=0.5)))
    return fig


def main():

    args = get_parser()
    result_dir = args.result_dir
    format = args.format

    with open(os.path.join(result_dir, 'y_test_true_list.pkl'), 'rb') as f, \
        open(os.path.join(result_dir, 'y_test_pred_list.pkl'), 'rb') as g:
        y_test_true_list = pickle.load(f)
        y_test_pred_list = pickle.load(g)
    y_true = [item for sublist in y_test_true_list for item in sublist]
    y_pred = [item for sublist in y_test_pred_list for item in sublist]

    performance = pd.read_csv(os.path.join(result_dir, 'performance.csv'))
    corrcoef_ave = float(performance['Corr_Coef'][0])
    corrcoef_std = float(performance['Corr_Coef'][1])

    mpl.style.use('seaborn-v0_8')
    fig = make_scatter_plot(y_true, y_pred, corrcoef_ave, corrcoef_std)
    outpath = f'{result_dir}/scatter.{format}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()
    print(f'Figure saved at {outpath}')


if __name__ == '__main__':
    main()