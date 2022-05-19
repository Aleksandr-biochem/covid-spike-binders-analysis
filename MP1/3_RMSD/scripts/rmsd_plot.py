import argparse
import os
import sys

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rcParams

def rmsd_graphs(trj, mutant, ax):
    """
    Plot RMSD ditributions of LCB, RBD and complex per each ns.

    param: trj - Object with RMSD metrics per trj (csv table)
    param: mutant - strain for which plot is generated (e.g. alpha+lcb1, omicron+lcb1)
    param: ax - layer of subplot

    """
    trj.columns = ['time_ns', 'Ca_lcb', 'Ca_rbd', 'Ca_complex']
    cols = list(trj.columns)
    n = len(cols)
    colors = ['slateblue', 'navy', 'red']
    for i in range(1, n):
        ax.plot(trj.time_ns, trj.iloc[:,i], label=cols[i], color=colors[i-1])
        
    ax.set_title(f'{mutant.capitalize()} strain RBD with LCB1', fontsize=45, fontweight="bold", pad=20)
    ax.set_xlabel('Time in ns', fontsize=35, fontweight="bold")
    #ax.set_xticks(np.arange(0, max(trj.time_ns) + 200, 200), fontsize=30)
    ax.set_ylabel(r'RMSD, $\bf{\AA}$', fontsize=35, fontweight="bold")
    #ax.set_yticks(np.arange(1.0, 3.25, 0.25), fontsize=30)
    
    ax.tick_params(axis='both', which='major', labelsize=35)

    ax.set_xlim(0, max(trj.time_ns))
    ax.set_ylim(0,3)
    ax.legend(fontsize=35, labels=[r"$sec.str. C\it{\alpha}\/\/LCB1\/\/\AA$",
        r"$sec. str. C\it{\alpha}\/\/RBD\/\/\AA$",
        r"$sec. str. C\it{\alpha}\/\/complex\/\/\AA$"],
        loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3)
    #plt.yticks(fontsize=30)
    #plt.xticks(fontsize=30)

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot rmsd distribution of LCB1, RBD and complex per each ns of simulated trajectories for all strains:')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing RMSD results', required=True, default="../../tables/RMSD")
    
    parser.add_argument('-o', nargs=1, help='Output filename', required=True)
    
    args = parser.parse_args()

    mutants = ['wt', 'alpha', 'delta', 'delta_plus', 'omicron']
    
    fig, ax = plt.subplots(len(mutants), 1, figsize = (65,40))
    ax = ax.flatten()

    for i, mutant in enumerate(mutants):
        data = pd.read_csv(os.path.join(args.i[0], f"{mutant}+mp1_rmsd_ca_ss.csv"))
        rmsd_graphs(data, mutant, ax[i])
    #plt.gca().set_color_cycle()
    #plt.yticks(fontsize=30)
    #plt.xticks(fontsize=30)
    plt.tight_layout()
    plt.savefig(f"../plots/{args.o[0]}", bbox_inches='tight')

    