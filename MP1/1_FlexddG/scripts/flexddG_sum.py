import os
import sys
from tqdm import tqdm
from glob import glob
import argparse

import pandas as pd
import numpy as np

def flexddG_data(rosetta_path, strain, score_func="fa_talaris2014"):
    """
    Summarize information of FlexddG data through generated subfolders.

    param: path - path to the dir with FlexddG results
    param: state - which state is summarized: wt-dG, mutant-dG, ddG
    param: backrup - number of backrups
    param: n_struct - number of structures to analyze
    param: score_func - which score is summarized
    
    """
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

    df_all = dict()
    rids_of_interest = range(1, 37)
    # Start data preparation
    for residue_number in tqdm(rids_of_interest):
        structures = (35, 35)
        path_to_saturation_csv = glob(f'{rosetta_path}/{residue_number:02d}_*/analysis_output_{strain}_pdb/output_saturation-results.csv')[0]
        data_ddG = pd.read_csv(path_to_saturation_csv)
        sign_col = ['score_type_name', 'backrub_steps', 'nstruct', 'case_name', 'score_function_name', 'scored_state', 'total_score']
        
        # extract ddG
        data_filt = data_ddG[data_ddG.scored_state.isin(['ddG'])][sign_col]
        
        # get only fa_talaris2014 function 
        data_filt_def_score = data_filt[(data_filt['score_function_name'] == score_func) & (data_filt.backrub_steps == 35000) & (data_filt.nstruct >= structures[0]) & (
                            data_filt.nstruct <= structures[1])]

        # sort data frame by total score and case name
        data_filtered_ds_sorted = data_filt_def_score.sort_values(["case_name", "total_score"], ascending=True)
        temp = [name.split("_")[-1] for name in data_filtered_ds_sorted["case_name"].to_numpy()]
        data_filtered_ds_sorted.insert(0, 'case_name_1', temp)
        ordered_classes = list(one_to_three_leter_name.keys())


        df_list = []

        for i in ordered_classes:
            df_list.append(data_filtered_ds_sorted[data_filtered_ds_sorted['case_name_1']==i])

        ordered_df = pd.concat(df_list)

        names = ordered_df.case_name_1.to_numpy()
        ddG = ordered_df["total_score"].to_numpy().astype(np.float64)
        #ordered_df.total_score = [float(i) for i in ordered_df.total_score]
        res = os.path.realpath(path_to_saturation_csv).split('/')[7]
        res = f"{res.split('_')[0]}_{three_to_one_leter_name[res.split('_')[1]]}"

        #res_df =  ordered_df.groupby('case_name_1', as_index = False).mean() 
        if 'residue_id' not in df_all:
            df_all['residue_id'] = names
        df_all[res] = ddG

    # Transform to dataframe format
    data = pd.DataFrame(df_all)
    data.set_index('residue_id', inplace = True)
    return data

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Summarize information about ddG from FlexddG output folders:')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing subfolders with FlexddG output', required=True, default="./MP1")

    parser.add_argument('-s', nargs=1,   help='Strain to analyze', required=True)
    
    parser.add_argument('-f', nargs=1,   help='Get values for score_function: fa_talaris2014 or fa_talaris2014-gam', required=True)
    
    parser.add_argument('-o', nargs=1,   help='Output .csv name', required=True)
    
    args = parser.parse_args()
    mutation_data = flexddG_data(rosetta_path=args.i[0], strain=args.s[0], score_func=args.f[0])
    mutation_data.to_csv(f'{args.o[0]}.csv')
