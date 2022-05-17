import sys
sys.path.append('../../../md_utils')
from md_utils.frame.extract import extract_one_letter_amino_acid_seq
from pyxmolpp2 import PdbFile, mName
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_ddg(df, path_to_pdb, fixed_mutation, mutant, backrub_steps, nstruct, score):

    # extract nonmutated sequence of MP3
    structure = PdbFile(path_to_pdb).frames()[0]
    MP_sequence = extract_one_letter_amino_acid_seq(structure, molname="A")

    #filter ddG data
    display_columns = ['backrub_steps', 'case_name', 'nstruct', 'score_function_name', 'scored_state', 'total_score', 'total_score.1']
    ddg_data = df.loc[ df['scored_state'] == 'ddG' ][display_columns]

    # filter data elements
    ddg_data_filt = ddg_data[(ddg_data['backrub_steps'] == backrub_steps) & (ddg_data['nstruct'] == nstruct) 
                             & (ddg_data['score_function_name'] == score)]

    # figures to plot
    mean_metrics = ddg_data_filt.loc[:, 'total_score'].astype(float).to_numpy().round(2)
    std_metrics = ddg_data_filt.loc[:, 'total_score'].astype(float).to_numpy()


    # a function to set axis parameters
    def set_axis_parameters(xname="", yname="", title_tag=""):
        fig = plt.figure(figsize=(14, 9))
        ax = fig.add_subplot(111)

        ax.set_xlabel('{xname}'.format(xname=xname), fontsize=25, weight="bold", style='italic', labelpad=15)
        ax.set_title('{title_tag}'.format(title_tag=title_tag), fontsize=35, loc='center', pad=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        return fig, ax

    # a function to label boxes
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            if height > 0:
                ax.text(rect.get_x() + rect.get_width()/1.6, height*1.2,
                        height,
                        ha='left', va='bottom', fontsize=15)
            else:
                ax.text(rect.get_x() + rect.get_width()/1.6, height*1.2,
                        height,
                        ha='left', va='bottom', fontsize=15)

    # plot
    fig, ax = set_axis_parameters()

    ind = np.arange(len(mean_metrics))  # the label locations
    width = 0.36  # the width of the bars

    rects = ax.bar(ind - width / 2, mean_metrics, width,
                   color="khaki", zorder=0, edgecolor="dimgrey", linewidth=2)

    ax.errorbar(ind - width / 2,
                mean_metrics,
                ls="",
                lw=2,
                yerr=std_metrics,
                capsize=5,
                color="black",
                zorder=5)

    # set ticks and labels
    ax.set_xticks(ind - width / 2)

    ax.set_xticklabels([f"{MP_sequence[int(case.split('_')[-1][:-1]) - 1]}{case.split('_')[-1]}" for case in ddg_data_filt['case_name']], fontsize=15)

    autolabel(rects)
    plt.axhline(y=0, color='dimgrey', alpha=0.6, linestyle='-', lw=2)
    ax.set_ylabel(f"{score} score", weight="bold", style='italic', fontsize=20, labelpad=10)

    ax.set_title(
        f"ΔΔG for pair mutations in MP3/{mutant} with {fixed_mutation}",
        fontsize=20, weight="bold", pad=15)

    ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize = 15)

    plt.savefig(f"../plots/{mutant}_pair_mutations.png", bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot pair mutations result from FlexddG scan:')
    
    parser.add_argument('-i', nargs=1, help='Path to input struct_scores_results.csv data', required=True)

    parser.add_argument('-p', nargs=1, help='Path to FlexddG input .pdb structure', required=True)

    parser.add_argument('-mutant', nargs=1, help='Specify RBD mutant name (alpha, delta+, etc.)', required=True) 

    parser.add_argument('-s', nargs=1, help='Fixed substitution, e.x. D37R', required=True) 

    parser.add_argument('-b', nargs=1, help='The backrup step number to summarize, integer', required=True) 

    parser.add_argument('-n', nargs=1,   help='Values for nstruct to filter, integer', required=True)    

    parser.add_argument('-f', nargs=1,   help='Value for score_function to filter: fa_talaris2014 or fa_talaris2014-gam', required=True)   
    
    args = parser.parse_args()

    # path to .csv data
    PATH_TO_DATA = args.i[0]
    # path to pdb
    path_to_pdb = args.p[0]
    # mutant analyzed
    mutant = args.mutant[0]
    # fixed mutation information
    fixed_mutation = args.s[0]
    # backrup steps state to summarize
    backrub_steps = int(args.b[0])
    # nstruct to filter
    nstruct = int(args.n[0])
    # score function
    score = args.f[0]

    aa_data = pd.read_csv(PATH_TO_DATA, index_col=0)
    plot_ddg(aa_data, path_to_pdb, fixed_mutation, mutant, backrub_steps, nstruct, score)