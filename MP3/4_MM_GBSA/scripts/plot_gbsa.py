import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def plot_mmgbsa(path_to_data, igb, start_frame, sate):

    # possible mutants for analysis
    mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

    # possible labels for complexes
    labels = {"wt+mp3": "Wild\ntype", "alpha+mp3": "Alpha",
              "delta+mp3": "Delta", "delta_plus+mp3": "Delta+", "omicron+mp3": "Omicron",
              "delta_p+mp3_d37r": "\nDelta+\nMP3(D37R)", "delta_p+mp3_d37r_t10w": "\n\n\nDelta+\nMP3(T10W;D37R)"}

    # lists to collect mean binding energies with std for igb2 and igb8 runs
    dG = []
    selected_labels = [] # labels for graphs

    for mutant in mutants:
        try:
            data = pd.read_csv(os.path.join(path_to_data, f"{mutant}_{igb}_salt150.csv"), index_col=0).iloc[start_frame:, :]
            dG.append([data['delta'].mean(), data['delta'].std()])
            selected_labels.append(labels[mutant])
        except:
            continue

    # function for boxplot labeling with mean energies
    def autolabel_ddG(rects):
        for rect in rects:
            height = rect.get_height()
            if height > 0:
                h = -1
            else: 
                h = 0.3
            plt.text(rect.get_x(), h,
                    height,
                    ha='left', va='bottom', fontsize=15)

    def autolabel_dG(rects):
        max_h = 0
        for rect in rects:
            height = rect.get_height()
            if height < max_h:
                max_h = height
            plt.text(rect.get_x(), height*1.35,
                height,
                ha='left', va='bottom', fontsize=15)
        plt.ylim(max_h-20, 0)

    # create boxplot
    fig = plt.figure(figsize=(10, 5))

    dG = np.array(dG)
    mean_metrics = dG[:, 0]
    std_metrics = dG[:, 1]

    # correct data according to plotting state
    if state == 'ddG':
        mean_metrics = np.array(mean_metrics)
        mean_metrics = mean_metrics[1:] - mean_metrics[0]
        selected_labels = selected_labels[1:]

    fig.add_subplot(1, 1, 1)
    ind = np.arange(len(mean_metrics))  # the label locations
    width = 0.36                        # the width of the bars

    rects = plt.bar(ind - width / 2, mean_metrics.round(2), width,
                    color="khaki", zorder=0, edgecolor="dimgrey", linewidth=2)

    if state == 'dG':
        plt.errorbar(ind - width / 2,
                    mean_metrics,
                    ls="",
                    lw=2,
                    yerr=std_metrics,
                    capsize=5,
                    color="black",
                    zorder=5)

    # set ticks and labels
    plt.xticks(ticks = (ind - width / 2), labels=selected_labels, fontsize=12)
    plt.yticks(fontsize=12)

    plt.axhline(y=0, color='dimgrey', alpha=0.6, linestyle='-', lw=2)

    if state == 'dG':
        autolabel_dG(rects)
        plt.ylabel(r"$\Delta G_{binding}$", fontsize=20)
        plt.title(f"MM-GBSA ΔG binding for RBD complexes with MP3 ({igb})", fontsize=15, pad=10)
    else:
        autolabel_ddG(rects)
        plt.ylabel(r"$\Delta\Delta G_{binding}$", fontsize=20)
        plt.title(f"MM-GBSA ΔΔG binding for RBD/MP3 complexes in relation to wild type ({igb})", fontsize=15, pad=10)

    # save image
    plt.savefig(f"../plots/{state}_mmgbsa_{igb}.png", bbox_inches='tight', dpi=300)

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot MM-GBSA results:')
    
    parser.add_argument('-i', nargs=1, help='Path to directory with .csv data', required=True) 
    
    parser.add_argument('-igb', nargs=1,   help='Specify igb run: igb2 or igb8', required=True)

    parser.add_argument('-f', nargs=1,   help='A frame to start data aggregation from, default 500', required=False)

    parser.add_argument('-s', nargs=1,   help='A state to plot: dG or ddG (computed versus wt complex)', required=False)
    
    args = parser.parse_args()

    # path to .csv data
    PATH_TO_DATA = args.i[0]
    # igb type to summarize
    igb = args.igb[0]
    # starting frame
    try:
        start_frame = int(args.f[0])
    except:
        start_frame = 500

    state = args.s[0]

    plot_mmgbsa(PATH_TO_DATA, igb, start_frame, state)