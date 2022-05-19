import os
import sys
import argparse

import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors
import seaborn as sns

# class for nirmalization of ddG values
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

# function to highlight cells
def highlight_cell(x, y, ax=None, **kwargs):
    rect = plt.Rectangle((x - .5, y - .5), 1, 1, fill=False, **kwargs)
    ax = ax or plt.gca()
    ax.add_patch(rect)
    return rect

# function to create heatmap
def heatmap_rosetta(summary, strain, output='.'):
    """
    Generate heatmap of primary LCB-to-mutant LCB dG difference.
    
    param: rosetta_path - path to the folder with FlexddG results
    param: strain - which strain is analyzed
    param: output - path to output folder

    """

    data_corrected = pd.read_csv(summary, index_col=0)
    # Assign subplot layers
    fig, ax = plt.subplots(figsize=(20, 10))

    vmin = np.array(data_corrected).min()
    vmax = np.array(data_corrected).max()
    norm = MidpointNormalize(vmin=vmin, vmax=vmax, midpoint=0)

    cmap = 'RdYlBu_r'
    img = plt.imshow(data_corrected, cmap=cmap, norm=norm)

    #col_labels = [col + '\n' + str(i)  else col.split('_')[1] for i, col in enumerate(data_corrected.columns) if i in range(1,36, 5)]
    col_labels = []
    scale_list = [j for j in range(0, 36, 5)]
    for i, col in enumerate(list(data_corrected.columns)):
        if i in scale_list:
            i2 = i + 1
            pat = col.split('_')[1] + '\n' + str(i2)
        else:
            pat = col.split('_')[1]
        col_labels.append(pat)

    #ax.set_xticks(np.arange(len(list(data_corrected.columns))), labels=list(data_corrected.columns), fontsize=15)
    ax.set_xticks(np.arange(len(list(data_corrected.columns))), labels=col_labels, fontsize=15)
    ax.set_yticks(np.arange(len(list(data_corrected.index))), labels=list(data_corrected.index), fontsize=15)

    ax.set_title(f'The saturation mutational scan of MP1 in complex with {strain} RBD', fontsize=20, pad=20, fontweight="bold")
    ax.set_xlabel('residues of MP1', fontsize=20, fontweight="bold")
    ax.set_ylabel('mutated residues of MP1', fontsize=20, fontweight="bold")

    #plt.setp(ax.get_xticklabels(), rotation = 45, ha="right",
    #         rotation_mode="anchor")

    # cbar settings
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.15)
    cbar = plt.colorbar(img, cax=cax)
    cbar.ax.tick_params(labelsize=20)

    # highlight cells
    value = -1
    for index_x, index_y in zip(np.where(data_corrected <= value)[-1], np.where(data_corrected <= value)[0]):
        highlight_cell(index_x, index_y, color="blue", linewidth=2.5, ax=ax)
        text = ax.text(index_x, index_y, round(data_corrected.iloc[index_y, index_x],1),
                           ha="center", va="center", color="black")

    # highlight overlapping amino acids
    correlations = np.where([[x in element for x in data_corrected.index.tolist()] for element in data_corrected.columns.tolist()])
    for index_y, index_x in zip(correlations[-1], correlations[0]):
        highlight_cell(index_x, index_y, color="grey", linewidth=2.5, ax=ax)
        text = ax.text(index_x, index_y, round(data_corrected.iloc[index_y, index_x], 1),
                       ha="center", va="center", color="black")


    plt.savefig(output)

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Plot heatmap for ddG of mutated complex against wt')
    
    parser.add_argument('--input', '-i', nargs=1, help='Path to rosetta analysis folder', required=True, default="./LCB1")
    
    parser.add_argument('--strain', '-s', nargs=1, help='Used name of strain', required=True, default="delta+")
    
    parser.add_argument('--output', '-o', nargs=1,   help='Output filename', required=True)
    
    args = parser.parse_args()
    heatmap_rosetta(summary=args.input[0], strain=args.strain[0], output=args.output[0])