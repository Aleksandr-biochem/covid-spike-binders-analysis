import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse

def rmsd_plot(trj, ax, title, what_to_plot):
    columns = list(trj.columns)[1:]
    plotting_columns = {1: [2, 1, 0], 2: [0, 3, 4]}
    
    colors = ['#1b1954', '#0400ff', '#f50000']
    alphas = [0.8, 0.5, 0.5]
    
    for j, c, a in zip(plotting_columns[what_to_plot], colors, alphas):
        ax.plot(trj.time_ns, trj[columns[j]], label=columns[j], color=c, alpha=a)

    ax.set_title(title, fontsize=25, pad=20)
    ax.set_xlabel('time, ns', fontsize=25, labelpad=10)
    ax.set_ylabel(r'RMSD, $\rm\AA$', fontsize=25, labelpad=10)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.set_ylim(1, 4)
    ax.set_xlim(0, len(trj.time_ns))
    ax.legend(fontsize=17)
    
# specify mutants for plotting
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

titles = ["Wild type RBD with MP3",
          "Alpha RBD with MP3",
          "Delta RBD with MP3",
          "Delta+ RBD with MP3",
          "Omicron RBD with MP3",
          "Delta+ RBD with MP3(D37R)",
          "Delta+ RBD with MP3(T10W;D37R)",]

plotting_states = {1: 'rbd_vs_mp', 2: 'mp_vs_helices'}

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot RMSD:')
    
    parser.add_argument('-i', nargs=1, help='Path to data folder', required=True)

    parser.add_argument('-v', nargs=1, help='1 - plot RMSD of complex vs RBD and MP; 2 - plot MP RMSD vs its helices', required=True)    
    
    args = parser.parse_args()

    path_to_data = args.i[0]
    what_to_plot = int(args.v[0])

    for title, mutant in zip(titles, mutants):
        try: 
            # read data 
            trj_data = pd.read_csv(os.path.join(path_to_data, f'rmsd_ca_ss_{mutant}.csv'))

            fig, ax = plt.subplots(1, 1, figsize=(20, 5))
            rmsd_plot(trj_data, ax, title, what_to_plot)
            plt.savefig(f"../plots/rmsd_{mutant}_{plotting_states[what_to_plot]}.png", bbox_inches='tight', dpi=300)
            plt.close(fig)
        except:
            continue

