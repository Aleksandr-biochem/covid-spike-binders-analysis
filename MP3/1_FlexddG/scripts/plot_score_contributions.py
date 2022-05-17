import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse

def plot_scores(data, old_residue, residue_id, new_residue, variant, backrub_steps):
    
    # a function to tune axis parameters
    def set_axis_parameters(xname="", yname="", title_tag=""):
        fig = plt.figure(figsize=(14, 9))
        ax = fig.add_subplot(111)

        ax.set_xlabel('{xname}'.format(xname=xname), fontsize=25, labelpad=15)
        ax.set_title('{title_tag}'.format(title_tag=title_tag), fontsize=35, loc='center', pad=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        return fig, ax

    # a function to label contribution values to bars
    def autolabel(rects):
        for rect in rects:
            height = rect.get_height()
            if height > 0:
                ax.text(rect.get_x() + rect.get_width()/1.6, height*1.1,
                        height,
                        ha='left', va='bottom', fontsize=15, fontweight="bold")
            else:
                ax.text(rect.get_x() + rect.get_width()/1.6, 0.1,
                        height,
                        ha='left', va='bottom', fontsize=15, fontweight="bold")
    
    # filter data from FlexddG output 
    case_name = data.iloc[0, 0][:-2]
    df_filtered = data[(data['case_name'] == f'{case_name}_{new_residue}') & 
                    (data['backrub_steps'] == backrub_steps)]
    
    # aggregate data on contributions
    metrics = ["fa_atr", "fa_elec", "fa_rep", "fa_sol", "hbond_bb_sc", "hbond_sc", "hbond_lr_bb"]
    mean_metrics = df_filtered[metrics].mean()
    std_metrics = df_filtered[metrics].std()

    # plot
    fig, ax = set_axis_parameters()

    ind = np.arange(len(metrics))  # the label locations
    width = 0.36  # the width of the bars
    
    # plot bars
    rects = ax.bar(ind - width / 2, mean_metrics.round(2), width,
                   color="khaki", zorder=0, edgecolor="dimgrey", linewidth=2)
    
    # plot errors
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
    ax.set_xticklabels([r"$\Delta$" + metric for metric in metrics], fontsize=15)

    autolabel(rects)
    plt.axhline(y=0, color='dimgrey', alpha=0.6, linestyle='-', lw=2)
    ax.set_ylabel("REU", fontsize=20)
    
    label = r"$\Delta\Delta G_{total}=$"
    ax.set_title(
        f"Energy terms MP3({old_residue}{residue_id}{new_residue}) with {variant}, {label}{mean_metrics.sum().round(2)}",
        fontsize=20, pad=15)

    # ax.set_ylim(-5, 5)
    
    plt.savefig(f"../plots/score_contributions_MP3_{old_residue}{residue_id}{new_residue}.png", bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot score contributions for a particular mutaion from FlexddG scan:')
    
    parser.add_argument('-i', nargs=1, help='Path to input struct_scores_results.csv data', required=True)

    parser.add_argument('-mutant', nargs=1, help='Specify RBD mutant name (alpha, delta+, etc.)', required=True) 

    parser.add_argument('-start_res', nargs=1, help='The residue of mutation in single letter code', required=True) 

    parser.add_argument('-pos', nargs=1, help='Position of mutation, integer', required=True) 

    parser.add_argument('-mut_res', nargs=1, help='The residue after mutation in single letter code', required=True) 

    parser.add_argument('-b', nargs=1, help='The backrup step number to summarize, integer', required=True)        
    
    args = parser.parse_args()

    # path to .csv data
    PATH_TO_DATA = args.i[0]
    # mutant analyzed
    mutant = args.mutant[0]
    # fixed mutation information
    initial_residue, index, mutated_residue = args.start_res[0], args.pos[0], args.mut_res[0] 
    # backrup steps state to summarize
    backrub_steps = int(args.b[0])

    aa_data = pd.read_csv(PATH_TO_DATA, index_col=0)
    plot_scores(aa_data, initial_residue, index, mutated_residue, mutant, backrub_steps)

