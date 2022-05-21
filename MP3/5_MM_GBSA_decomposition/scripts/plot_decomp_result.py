import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse
from os import path

params = {
    'font.size': 42,
    'axes.labelsize': 70,
    'axes.labelpad': 0,
    'axes.titlepad': 140,
    'axes.titleweight': 'bold',
    'axes.edgecolor': 'black',
    'axes.titlesize': 64,
    'figure.max_open_warning': 2,
    'figure.figsize': (60, 40),
    'legend.fontsize': 40,
    'legend.handlelength': 2,
    'xtick.labelsize': 38,
    'ytick.labelsize': 38,
    'axes.axisbelow': True
}

SUBPLOTTITLEFONTSIZE = 50
MARGIN = 0.02
BARWIDTH = 0.8
plt.rcParams.update(params)


def heatmap(data, row_labels, col_labels, ax=None, cbarlabel="", **kwargs):
    if not ax:
        ax = plt.gca()

    # generate the heatmap
    im = ax.imshow(data, **kwargs)

    # create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", fontsize=55, labelpad=20)

    # set ticks and label
    n_rows, n_cols = data.shape
    ax.set_yticks(np.arange(n_rows), labels=row_labels)
    ax.set_xticks(np.arange(n_cols), labels=col_labels)

    # set non-visible lines confining the plot area
    ax.spines[:].set_visible(False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.1f}",
                     textcolors=("black", "white"),
                     threshold=None, vmin=-32, vmax=32, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    im.set_clim(vmin=vmin, vmax=vmax)

    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = np.abs(vmin) / 2. - 0.5

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if np.abs(data[i, j]) >= 1:
                kw.update(color=textcolors[int(np.abs(data[i, j]) > threshold or np.abs(data[i, j]) < 1)],
                          fontsize=36)
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

    return texts


def set_axis_parameters(xname="", yname="", title_tag=""):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('{xname}'.format(xname=xname))
    ax.set_ylabel('{yname}'.format(yname=yname))
    ax.set_title('{title_tag}'.format(title_tag=title_tag), loc='center')
    return fig, ax


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, rId
    import pandas as pd
    import os

    # parse arguments
    parser = argparse.ArgumentParser(description='Plot MM-GBSA decomposition results as a differential contact maps:')
    
    parser.add_argument('-i', nargs=1, help='Path to directory with .csv data', required=True) 
    
    parser.add_argument('-igb', nargs=1,   help='Specify igb run to summarize: igb2 or igb8', required=True)
    
    args = parser.parse_args()

    path_to_data = args.i[0]
    igb = args.igb[0]

    # specify mutants for the analysis
    mutants = ["alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

    ### read wt structure
    path_to_wt = f"../../2_MD_Amber/wt+mp3/0_prepare/protein_named.pdb"
    ref_wt = PdbFile(path_to_wt).frames()[0]

    # read decomposition data for wt and initiate union residue parirs set
    df_decomp_wt = pd.read_csv(os.path.join(path_handling_decomp_dir, f"wt+mp3_{igb}_salt150_decomp.csv"))
    residue_pairs_wt = list(zip(df_decomp_wt["rId_mp"].values, df_decomp_wt["rId_rbd"].values))
    union_residue_pairs = set(residue_pairs_wt) 

    ## gather union, overlap and unique residues for all complexes
    print("Analyzing interacting residues")
    overlapping_desidues = dict() # collect overlaping residues sets for each mutant-wt pair
    found_mutants = []
    for mutant in mutants:
        # check  if data for mutant exists
        if not path.exists(os.path.join(path_to_data, f"{mutant}_{igb}_salt150_decomp.csv")):
            continue
        else: 
            found_mutants.append(mutant)  

            # read parsed decomp data with residue pairs, whose trajectory-average energy exceeded 1 kcal/mol
            df_decomp_mutant = pd.read_csv(os.path.join(path_handling_decomp_dir, f"decomp_{mutant}_{igb}_salt150_decomp.csv"))
            residue_pairs_mutant = list(zip(df_decomp_mutant["rId_mp"].values, df_decomp_mutant["rId_rbd"].values))

            # udate union residues in wt and mutant complexes and extract overlap for this exact pair
            union_residue_pairs = union_residue_pairs | set(residue_pairs_mutant)
            overlapping_desidues[mutant] = list(set(residue_pairs_wt) & set(residue_pairs_mutant))

    # process residues in union set 
    union_residue_mp, union_residue_rbd = map(list, zip(*union_residue_pairs))
    union_residue_mp, union_residue_rbd = np.unique(union_residue_mp), np.unique(union_residue_rbd)
    union_residue_mp.sort()
    union_residue_rbd.sort()

    # plot differenctial contact map for each mutant
    for mutant in found_mutants:
        print(f"Plotting {mutant}, {igb}")

        # read mutant reference
        path_to_mutant = f"../../2_MD_Amber/{mutant}/0_prepare/protein_named.pdb"
        ref_mutant = PdbFile(path_to_mutant).frames()[0]

        ### read parsed decomp data with residue pairs, whose trajectory-average energy exceeded 1 kcal/mol
        df_decomp_mutant = pd.read_csv(os.path.join(path_handling_decomp_dir, f"decomp_{mutant}_{igb}_salt150_decomp.csv"))
        residue_pairs_mutant = list(zip(df_decomp_mutant["rId_mp"].values, df_decomp_mutant["rId_rbd"].values))

        # get unique residues and overlap for this particular pair of complexes
        overlap_residue_pairs = overlapping_desidues[mutant]
        unique_residue_pairs_mutant = list(set(residue_pairs_mutant) - set(overlap_residue_pairs))
        unique_residue_pairs_wt = list(set(residue_pairs_wt) - set(overlap_residue_pairs))

        ### generate data_matrix for plot heatmap
        data_matrix = np.zeros((len(union_residue_mp), len(union_residue_rbd)))
        for ind_1, rid_mp in enumerate(union_residue_mp):
            for ind_2, rid_rbd in enumerate(union_residue_rbd):
                if (rid_mp, rid_rbd) in overlap_residue_pairs:
                    total_energy_wt = df_decomp_wt["total_energy_avg"][
                        (df_decomp_wt["rId_mp"] == rid_mp) & (df_decomp_wt["rId_rbd"] == rid_rbd)].values[0]
                    total_energy_mutant = df_decomp_mutant["total_energy_avg"][
                        (df_decomp_mutant["rId_mp"] == rid_mp) & (df_decomp_mutant["rId_rbd"] == rid_rbd)].values[0]
                    data_matrix[ind_1, ind_2] = total_energy_mutant - total_energy_wt
                elif (rid_mp, rid_rbd) in unique_residue_pairs_mutant:
                    total_energy_mutant = df_decomp_mutant["total_energy_avg"][
                        (df_decomp_mutant["rId_mp"] == rid_mp) & (df_decomp_mutant["rId_rbd"] == rid_rbd)].values[0]
                    data_matrix[ind_1, ind_2] = total_energy_mutant
                elif (rid_mp, rid_rbd) in unique_residue_pairs_wt:
                    total_energy_wt = df_decomp_wt["total_energy_avg"][
                        (df_decomp_wt["rId_mp"] == rid_mp) & (df_decomp_wt["rId_rbd"] == rid_rbd)].values[0]
                    data_matrix[ind_1, ind_2] = - total_energy_wt
                else:
                    data_matrix[ind_1, ind_2] = 0 # if residue_pair doesn't exist in df_decomp the corresponding total_eneregy is asssigned zero

        ### define rnames wt and mutant for labelling heatmap axes
        rnames_wt_rbd = [res.name for res in ref_wt.residues.filter(rId.is_in(set(union_residue_rbd)))]
        rnames_mutant_rbd = [res.name for res in ref_mutant.residues.filter(rId.is_in(set(union_residue_rbd)))]

        rnames_wt_mp = [res.name for res in ref_wt.residues.filter(rId.is_in(set(union_residue_mp)))]
        rnames_mutant_mp = [res.name for res in ref_mutant.residues.filter(rId.is_in(set(union_residue_mp)))]

        col_labels = [f'{r_name1}\n{r_num}\n{r_name2}' if r_name1 != r_name2 else f'{r_name1}\n{r_num}'
                      for r_name1, r_num, r_name2 in
                      zip(rnames_wt_rbd,
                          union_residue_rbd,
                          rnames_mutant_rbd)]

        row_labels = [f'{r_name1}{r_num}{r_name2}' if r_name1 != r_name2 else f'{r_name1}{r_num}'
                      for r_name1, r_num, r_name2 in
                      zip(rnames_wt_mp,
                          union_residue_mp,
                          rnames_mutant_mp)]

        # check HIS-HIE substitution
        for i in range(len(row_labels)):
            if 'HIS' in row_labels[i] and 'HIE' in row_labels[i]:
                row_labels[i] = row_labels[i][:-3]

        # check CYS-CIX substitution
        for i in range(len(col_labels)):
            if 'CYS' in col_labels[i] and 'CYX' in col_labels[i]:
                col_labels[i] = col_labels[i][:-3]

        ### define the common heatmap settings
        colormap = 'bwr'
        vmax = 15
        xname = '\nRBD-S'
        ligand = mutant.split('+')[1]
        variant = mutant.split('+')[0]
        if variant == 'delta_p' or variant == 'delta_plus':
            variant = 'delta+'
        ligand_for_title = ligand.upper()
        title_tag = f"Differential contact map: {variant}/{ligand_for_title} vs. wt/MP3 ({igb})"

        ### plot and save heatmap
        fig, ax = set_axis_parameters(xname=xname, yname=f"{ligand_for_title}\n", title_tag=title_tag)

        im, cbar = heatmap(data_matrix, row_labels=row_labels, col_labels=col_labels, ax=ax, cmap=colormap,
                           cbarlabel="ENERGY\n\n")

        annotate_heatmap(im, vmax=32)

        plt.savefig(f"../contact_map_{mutant}_{igb}.png", bbox_inches='tight', dpi=300)
        plt.close(fig)
