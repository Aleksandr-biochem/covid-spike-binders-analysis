import os
import sys
import argparse

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from MMPBSA_mods import API as MMPBSA_API
from MMPBSA_mods import amber_outputs as amber


def polar_graphs(data, mutant="", res_pair="", outpth="."):
    cols = list(data.columns)
    n = len(cols)
    
    for i in range(1, n):
        plt.plot(data.time, data.iloc[:,i], label=cols[i])
    res_pair_title = res_pair.upper()    
    plt.title(f"Polar interactions contribution into {res_pair_title} total binding energy of MP1 and {mutant.split('+')[0]} RBD-S", fontsize=37, fontweight="bold", pad=30)
    plt.xlabel('Time in ns', fontsize=35, fontweight="bold")
    plt.xticks(fontsize=30)
    plt.ylabel(r'$\bf{\Delta G_{binding}}$', fontsize=35, fontweight="bold")
    plt.yticks(fontsize=30)
    plt.xlim(500, data.time.max())
    #plt.set_ylim(0,3)
    #plt.legend(loc="upper left", bbox_to_anchor=(1,1), fontsize=20, labels=["Total binding \nenergy", "Polar interactions \nbinding energy"])
    plt.legend(fontsize=30, labels=["Total binding \nenergy", "Polar interactions \nbinding energy"])
    plt.tight_layout()
    plt.savefig(os.path.join(outpth, f"{mutant}_{res_pair}.png"))
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Plot polar interaction contribution into total dG per trj')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing full MM-GBSA decomposition data for residue of interest in .tsv format', required=True, default="../tables/")

    parser.add_argument('-r', nargs=1, help='Desired pair of residues to plot in the format AA№-AA№ where AA - amino acid in three letter format', required=True)

    parser.add_argument('-o', nargs=1, help='Path to output folder', required=False, default="./../plots/")
    
    args = parser.parse_args()

    mutants = ["alpha", "delta", "delta_plus", "omicron"]
    # residues = ["ASP30-ARG403", "ASP30-ASN417", "SER29-ASN417"]
    pair = args.r[0]
    
    for mut in mutants:
        plt.rcParams["figure.figsize"] = (35,15)    
        data = pd.read_csv(os.path.join(args.i[0],f"{mut}+mp1_polar_{pair}.tsv"), sep='\t')
        polar_graphs(data, mutant=mut, res_pair=pair, outpth=args.o[0])
        plt.clf()