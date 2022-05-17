import pandas as pd
from matplotlib.figure import Figure
import matplotlib as mpl
import numpy as np


three_to_one_leter_name = {
                    "GLY": "G",
                    "ALA": "A",
                    "VAL": "V",
                    "LEU": "L",
                    "ILE": "I",
                    "SER": "S",
                    "THR": "T",
                    "MET": "M",
                    "CYS": "C",
                    "ASN": "N",
                    "GLN": "Q",
                    "LYS": "K",
                    "ARG": "R",
                    "HIS": "H",
                    "ASP": "D",
                    "GLU": "E",
                    "PHE": "F",
                    "PRO": "P",
                    "TRP": "W",
                    "TYR": "Y",
                }

one_to_three_leter_name = {
                    "G": "GLY",
                    "A": "ALA",
                    "V": "VAL",
                    "L": "LEU",
                    "I": "ILE",
                    "S": "SER",
                    "T": "THR",
                    "M": "MET",
                    "C": "CYS",
                    "N": "ASN",
                    "Q": "GLN",
                    "K": "LYS",
                    "R": "ARG",
                    "H": "HIS",
                    "D": "ASP",
                    "E": "GLU",
                    "P": "PRO",
                    "F": "PHE",
                    "W": "TRP",
                    "Y": "TYR",
                }

structures = (35, 35)
initial_residue = "SER"
residue_number = 29
var = "pdb"

df = pd.read_csv('./analysis_output_deltaplus_pdb/output_saturation-results.csv')

df_filtered = df[
    (df.backrub_steps == 35000) & (df.score_function_name == "fa_talaris2014") & (df.nstruct >= structures[0]) & (
                df.nstruct <= structures[1]) & (df.scored_state == "ddG")]

df_filtered_sorted = df_filtered.sort_values(["case_name", "total_score", "total_score.1"], ascending=True)
temp = [name.split("_")[-1] for name in df_filtered_sorted["case_name"].to_numpy()]
df_filtered_sorted.insert(0, 'case_name_1', temp)
ordered_classes = list(one_to_three_leter_name.keys())

df_list = []

for i in ordered_classes:
   df_list.append(df_filtered_sorted[df_filtered_sorted['case_name_1']==i])

ordered_df = pd.concat(df_list)

names = ordered_df.case_name_1.to_numpy()
ddG = ordered_df["total_score"].to_numpy().astype(np.float64)
ddG_std = ordered_df[("total_score.1")].to_numpy().astype(np.float64)
names = list(map(lambda x: "{}".format(one_to_three_leter_name[x[-1]]), names))


df_filtered_zemu = df[
    (df.backrub_steps == 35000) & (df.score_function_name == "fa_talaris2014-gam") & (df.nstruct >= structures[0]) & (
                df.nstruct <= structures[1]) & (df.scored_state == "ddG")]

df_filtered_zemu_sorted = df_filtered_zemu.sort_values(["case_name", "total_score", "total_score.1"], ascending=True)
temp = [name.split("_")[-1] for name in df_filtered_zemu_sorted["case_name"].to_numpy()]
df_filtered_zemu_sorted.insert(0, 'case_name_1', temp)
ordered_classes = list(one_to_three_leter_name.keys())

df_list = []

for i in ordered_classes:
   df_list.append(df_filtered_zemu_sorted[df_filtered_zemu_sorted['case_name_1']==i])

ordered_df_zemu = pd.concat(df_list)

ddG_zemu = ordered_df_zemu["total_score"].to_numpy().astype(np.float64)
ddG_zemu_std = ordered_df_zemu[("total_score.1")].to_numpy().astype(np.float64)


fig = Figure(linewidth=4, figsize=(8, 4), dpi=300)
gs = mpl.gridspec.GridSpec(nrows=1, ncols=1, left=0.07, right=0.82, bottom=0.13, top=0.90, wspace=0.2, hspace=0.2)
ax1 = fig.add_subplot(gs[0, 0])

bar_width = 0.3
space_koeff = 0.5
ax1.bar(list(map(lambda n: n - space_koeff * bar_width, range(len(names)))), ddG, bar_width,
        yerr=ddG_std / structures[0] ** space_koeff, color='#92000a', capsize=3, label='ddG')
ax1.bar(list(map(lambda n: n + space_koeff * bar_width, range(len(names)))), ddG_zemu, bar_width,
        yerr=ddG_zemu_std / structures[0] ** space_koeff, color='#00008b', capsize=3, label='ddG_zemu')

ax1.set_xticks(list(map(lambda n: n, range(len(names)))))
ax1.set_xticklabels(names)
for tick in ax1.get_xticklabels():
    tick.set_rotation(60)
    tick.set_fontsize(10)
    tick.set_horizontalalignment('center')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
ax1.set_ylabel(u'\u0394\u0394G, REU', fontsize=12)
ax1.axhline(0, color='grey', linewidth=0.8)
ax1.set_title("{}{}".format(initial_residue, residue_number), fontsize=12, pad=6)
ax1.set_ylim(-1.5, 4.5)
fig.savefig("{}{}.png".format(initial_residue, residue_number),bbox_inches='tight')
