import argparse
import os
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


def get_parser():
    parser = argparse.ArgumentParser(
        usage=f'ex.) python {os.path.basename(__file__)}'
    )
    parser.add_argument(
        '--data', type=str, required=True,
        )
    parser.add_argument(
        '--outdir', type=str, required=True,
        )
    parser.add_argument(
        '--format', type=str, default='png',
        )
    return parser.parse_args()


def make_histgram(values):
    axis_interval = math.ceil((max(values)-min(values))/10)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(values, bins=50, color='gray', alpha=0.7)
    ax.set_xlabel('Value', fontsize=20, labelpad=15, weight='normal')
    ax.set_ylabel('Count', fontsize=20, labelpad=15, weight='normal')
    ax.xaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.yaxis.set_tick_params(direction='out', labelsize=15, width=3, pad=10)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(axis_interval))
    # ax.set_xlim(-25, 225)
    # ax.hist(values, color='gray', alpha=0.7, bins=1000)
    ax.text(0.8, 0.8, f"n={len(values)}", transform=ax.transAxes,
            ha='center', va='center', fontsize=15, alpha=1,
            bbox=(dict(boxstyle='square', fc='white', ec='gray', lw=0.5)))
    return fig


def main():

    args = get_parser()
    data_path = args.data
    result_dir = args.outdir
    format = args.format

    data = pd.read_csv(data_path)
    values = data['value'].values

    mpl.style.use('seaborn-v0_8')
    fig = make_histgram(values)
    outpath = f'{result_dir}/data_distribution.{format}'
    fig.savefig(outpath, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()
    print(f'Figure saved at {outpath}')


if __name__ == '__main__':
    main()
