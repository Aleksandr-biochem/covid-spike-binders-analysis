import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import os
import sys
import argparse

params = {
    'font.size': 42,
    'axes.labelsize': 48,
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
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

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
    ax.set_xlabel('{xname}'.format(xname=xname), fontweight="bold")
    ax.set_ylabel('{yname}'.format(yname=yname), fontweight="bold")
    ax.set_title('{title_tag}'.format(title_tag=title_tag), loc='center', fontweight="bold")
    return fig, ax


if __name__ == "__main__":
    from pyxmolpp2 import PdbFile, rId
    import pandas as pd
    import os
    
    parser = argparse.ArgumentParser(description='Plot difference matrix of mutant and wt dG ofinteraction on interface')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing MM-GBSA decomposition summary tables', required=True, default="../tables/")
    
    parser.add_argument('-w', nargs=1, help='Path to wt prepared structure with reasigned residue IDs in .pdb format', required=True)

    parser.add_argument('-t', nargs=1, help='Path to trj analysis folder', required=True)

    parser.add_argument('-o', nargs=1, help='Path to output folder', required=False, default="../plots/")
    
    args = parser.parse_args()
    
    ### read wt and mutant references
    #path_to_wt = f"/home/xenia/cov2/trj/mp1/wt+mp1_2/0_prepare/protein_named.pdb"
    path_to_wt = args.w[0]
    ref_wt = PdbFile(path_to_wt).frames()[0]

    mutants = ["alpha", "delta", "delta_plus", "omicron"]

    for mut in mutants:
        if mut == "omicron":
            path_to_mutant = os.path.join(args.t[0],f"{mut}+mp1_2/0_prepare/protein_named.pdb")
        else:
            path_to_mutant = os.path.join(args.t[0], f"{mut}+mp1/0_prepare/protein_named.pdb")
        ref_mutant = PdbFile(path_to_mutant).frames()[0]

        ### read parsed decomp data with residue pairs, whose trajectory-average energy exceeded 1 kcal/mol
        #path_decomp_dir = "/home/xenia/cov2/trj/MP1/summary_gbsa_decomp"
        path_decomp_dir = args.i[0]
        df_decomp_wt = pd.read_csv(os.path.join(path_decomp_dir, "wt+mp1_igb8_sum_decomp_GBSA.TDC.significant.tsv"), sep='\t', index_col=False)
        df_decomp_mutant = pd.read_csv(os.path.join(path_decomp_dir, f"{mut}+mp1_igb8_sum_decomp_GBSA.TDC.significant.tsv"), sep='\t', index_col=False)

        #df_decomp_wt = pd.read_csv(os.path.join(path_decomp_dir, "wt+mp1_igb8_sum_decomp_GBSA.TDC.significant.tsv"), sep='\t')
        #df_decomp_mutant = pd.read_csv(os.path.join(path_decomp_dir, "omicron+mp1_igb8_sum_decomp_GBSA.TDC.significant.tsv"), sep='\t')


        ### exctract union, overlap and unique residues in wt and mutant complex. It needs for generate data_matrix
        residue_pairs_wt = list(zip(df_decomp_wt["rId_mp"].values, df_decomp_wt["rId_rbd"].values))
        residue_pairs_mutant = list(zip(df_decomp_mutant["rId_mp"].values, df_decomp_mutant["rId_rbd"].values))

        union_residue_pairs = set(residue_pairs_wt) | set(residue_pairs_mutant)
        union_residue_mp, union_residue_rbd = map(list, zip(*union_residue_pairs))
        union_residue_mp, union_residue_rbd = np.unique(union_residue_mp), np.unique(union_residue_rbd)
        union_residue_mp.sort()
        union_residue_rbd.sort()

        overlap_residue_pairs = list(set(residue_pairs_wt) & set(residue_pairs_mutant))

        unique_residue_pairs_wt = list(set(residue_pairs_wt) - set(overlap_residue_pairs))
        unique_residue_pairs_mutant = list(set(residue_pairs_mutant) - set(overlap_residue_pairs))

        ### generate data_matrix for plot heatmap
        ### if residue_pair doesn't exist in df_decomp the corresponding total_eneregy is asssigned zero
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

        ### define rnames wt and mutant for labelling heatmap axes
        rnames_wt_rbd = [res.name for res in ref_wt.residues.filter(rId.is_in(set(union_residue_rbd)))]
        rnames_mutant_rbd = [res.name for res in ref_mutant.residues.filter(rId.is_in(set(union_residue_rbd)))]

        rnames_mutant_mp = [res.name for res in ref_wt.residues.filter(rId.is_in(set(union_residue_mp)))]

        col_labels = [f'{r_name1}\n{r_num}\n{r_name2}' if r_name1 != r_name2 else f'{r_name1}\n{r_num}'
                      for r_name1, r_num, r_name2 in
                      zip(rnames_wt_rbd,
                          union_residue_rbd,
                          rnames_mutant_rbd)]

        row_labels = [f"{n}{k}" for k, n in
                      zip(union_residue_mp,
                          rnames_mutant_mp)]

        ### define the common heatmap settings
        colormap = 'bwr'
        vmax = 15
        xname = '\nRBD-S'
        ligand = "mp1"
        ligand_for_title = ligand.upper()
        mutants_dict_title = {
        "alpha+mp1": "alpha",
        "delta+mp1": "delta",
        "delta_plus+mp1": r"$\bf{delta^+}$",
        "omicron+mp1": "omc"
        }
        title_tag = f"DIFFERENTIAL CONTACT MAP: {mutants_dict_title[f'{mut}+mp1']}-S / {ligand_for_title} vs. wt-S / {ligand_for_title}\n"

        ### plot and save heatmap
        fig, ax = set_axis_parameters(xname=xname, yname=f"{ligand_for_title}\n", title_tag=title_tag)

        im, cbar = heatmap(data_matrix, row_labels=row_labels, col_labels=col_labels, ax=ax, cmap=colormap,
                           cbarlabel="ENERGY\n\n")
        annotate_heatmap(im, vmax=32)

        out_dir = args.o[0]
        plt.savefig(os.path.join(out_dir, f"{mut.split('+')[0]}_vs_wt_igb8.png"))
