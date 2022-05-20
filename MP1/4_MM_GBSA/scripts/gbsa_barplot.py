import os
import sys
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from MMPBSA_mods import API as MMPBSA_API
from MMPBSA_mods import amber_outputs as amber

def set_axis_parameters(df, fig=None, ax=None, run=""):
    """
    Set axis parameters and labels.

    param: df - basic frame
    param: fig - fig object of announced plot
    param: ax - axis object of announced plot
    param: run - number of solvate model
    """
    if (fig is not None and ax is not None):
        x_ticks = ['Alpha', 'Delta', r'Delta$^+$', 'Omicron']
        ax.set_xticks(range(0,4), x_ticks)
        ax.tick_params(axis='both', which='major', labelsize=20)        
        ax.set_ylabel(r"$\bf{\Delta \Delta G}$", fontsize=20, fontweight="bold")
        ax.set_title("Difference of free binding energy (\u0394G) \n between wt and mutant MP1-RBD complexes under igb{0} model".format(run), fontweight="bold", fontsize=20, pad=30)
        return fig, ax
    else:
        print("Plot was not assigned for attaching subplots. Assign it first.")

def autolabel(ax=None, rects=None, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0, 'right': 1, 'left': -1}

    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(round(height, 2)),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(offset[xpos]*3, 3),  # use 3 points offset
                    textcoords="offset points",  # in both directions
                    ha=ha[xpos], va='bottom', fontsize=20)

def plot_gbsa_bar(path=".", outpth=".", run=""):
    """
    Plot barplots of difference between wt complex and mutant complex MM-GBSA.

    param: path - path to directory with MM-GBSA summary tables
    param: outpth - path to output directory
    """

    complexes = ['alpha', 'delta', 'delta_plus', 'omicron']
    #gbsa_path = "/home/xenia/cov2/trj/LCB1/summary_gbsa/"
    gbsa_path = path
    fig, ax = plt.subplots(1,figsize=(15, 10))
    # prepare dataframe
    df = pd.DataFrame(columns=complexes)
    # load data from all complexes per igb model
    data_wt = pd.read_csv(os.path.join(gbsa_path, f"wt+mp1_igb{run}_salt150.tsv"), sep="\t", index_col=0)
    for cmplx in complexes:
        data = pd.read_csv(os.path.join(gbsa_path, f"{cmplx}+mp1_igb{run}_salt150.tsv"), sep="\t", index_col=0)
        df[cmplx] = abs(data_wt.dG[500:]) - abs(data.dG[500:])
    # obtain mean and std for bars' vis.
    value = df.mean()
    std = df.std()
    # start subplotting
    ax.axhline(y=0, zorder=0)
    # draw bars 
    rects = ax.bar(np.arange(len(df.columns)), value, width=0.2, 
        align='center', alpha=0.5, capsize=10,
        color = "khaki", edgecolor="dimgrey", linewidth=2)
    fig, ax = set_axis_parameters(df, fig, ax, run=run)
    autolabel(ax=ax, rects=rects, xpos="right")
            
    plt.tight_layout()
    plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)
    plt.savefig(os.path.join(outpth, f"MP1_dG_GBSA_igb{run}.png"))

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot difference of wt and mutant complex MM-GBSA:')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing summary statistics of MM-GBSA', required=True, default="../tables/")

    parser.add_argument('-r', nargs=1, help='Number of model to analyze', type=str, required=True, default="2")
    
    parser.add_argument('-o', nargs=1, help='Path to output folder', required=True, default="../plots")
    
    args = parser.parse_args()

    plot_gbsa_bar(path=args.i[0], outpth=args.o[0], run=args.r[0])
