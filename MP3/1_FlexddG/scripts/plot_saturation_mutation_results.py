import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rcParams
import argparse

def plot_heatmap(PATH_TO_DATA, mutant, output='.'):

    # load data 
    data_deltap = pd.read_csv(PATH_TO_DATA, index_col=0)
    data_deltap = data_deltap.sort_values(by='id')
    data_deltap = data_deltap.reset_index(drop=True)

    # class to adjust color map
    class MidpointNormalize(colors.Normalize):
        def __init__(self, vmin, vmax, midpoint=0, clip=False):
            self.midpoint = midpoint
            print(vmin, vmax)
            colors.Normalize.__init__(self, vmin, vmax, clip)

        def __call__(self, value, clip=None):
            normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
            normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
            normalized_mid = 0.5
            x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
            return np.ma.masked_array(np.interp(value, x, y))

    # set up figure
    plt.figure(figsize=(40, 17))

    # get ddG data to plot
    cols = list(data_deltap.columns)[2:] # columns to plot
    ddg = data_deltap[cols].to_numpy().transpose()

    # annotate only cells with values below -1, these mutations are considered stabilyzing 
    an = np.vectorize(lambda x: ' ' if x>=-1 else str(round(x, 1)))(ddg)

    # set heatmap parameters
    vmin, vmax = -3, 4 # min and max values on color bar
    norm = MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0)
    cmap = 'RdYlBu_r'

    # write x labels
    x_labels = []
    for i, r in zip(data_deltap.id, data_deltap.residue):
        if i == 1 or i % 5 == 0:
            x_labels.append(f'{r}\n{i}')
        else:
            x_labels.append(f'{r}')

    # set up heatmap
    s = sns.heatmap(data=ddg, 
                    xticklabels=x_labels,
                    yticklabels=cols, 
                    vmin=vmin, vmax=vmax, center=0,
                    annot=an,
                    annot_kws = {'size':20, 'rotation': 45, "weight": "bold"},
                    fmt = '',
                    norm=norm,
                    cmap=cmap,
                    cbar_kws={"pad": 0.01})

    # highlight substitutions to the same amino acid in drey dashed patches
    residues = list(data_deltap.residue)
    for i in range(len(residues)):
        for j in range(len(cols)):
            if residues[i] == cols[j]:
                s.add_patch(plt.Rectangle((i, j), 1, 1, fill=False,
                            linestyle='--', alpha=0.5,
                            edgecolor='grey', lw=5))
                
    # highlight stabilizing substitutions in green patches
    for i in range(ddg.shape[0]):
        for j in range(ddg.shape[1]):
            if ddg[i, j] <= -1:
                s.add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                             edgecolor='#2dcd00', lw=5))
              
                
    # set labels
    cbar = s.collections[0].colorbar
    cbar.ax.tick_params(labelsize=35)

    plt.title(f"FlexddG ΔΔG scores for MP3 with {mutant} variant in relation to wild type", weight="bold", fontsize=45, pad=40)

    s.set_xlabel('MP3 residues', fontsize=40, style='italic', labelpad=20)
    s.set_ylabel('Residue after mutation', style='italic', fontsize=40, labelpad=20)

    s.set_yticklabels(s.get_ymajorticklabels(), fontsize = 30, rotation=0)
    s.set_xticklabels(s.get_xmajorticklabels(), fontsize = 30)


    plt.savefig(f'{output}/heatmap_{mutant}.png', bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot heatmap of FlexddG saturation mutational scan:')
    
    parser.add_argument('-i', nargs=1, help='Path to input .csv data', required=True)

    parser.add_argument('-mutant', nargs=1, help='Specify RBD mutant name (alpha, delta+, etc.)', required=True)    
    
    parser.add_argument('-o', nargs=1,   help='Path to output figure.png', required=True)
    
    args = parser.parse_args()

    # path to .csv data
    PATH_TO_DATA = args.i[0]
    # mutant analyzed
    mutant = args.mutant[0]
    # otput dir
    output = args.o[0]

    plot_heatmap(PATH_TO_DATA, mutant, output)
