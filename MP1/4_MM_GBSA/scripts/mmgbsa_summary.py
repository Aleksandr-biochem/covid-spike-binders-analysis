from MMPBSA_mods import API as MMPBSA_API
import pickle

import pandas as pd
import numpy as np
import os
import argparse
# the workflow saves data in csv for subsequent usage

def pickle_mmgbsa_data(path_to_info_file, path_to_output_dir="./", cmplx=None, run=None):
    wd = os.getcwd()
    os.chdir(os.path.dirname(path_to_info_file))
    data = MMPBSA_API.load_mmpbsa_info("_MMPBSA_info")
    os.chdir(wd)
    os.makedirs(path_to_output_dir, exist_ok=True)
    with open(os.path.join(path_to_output_dir, f'MMPBSA_data_{cmplx}+mp1_{run}.pickle'), 'wb') as f:
        pickle.dump(data, f)

def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
    with open(pcike_file, 'rb') as f:
        data = pickle.load(f)
    return data

def sum_gbsa(input_dir=".", out_dir="summary_gbsa", answ="y"):
     complexes = ["wt", "alpha", "delta", "delta_plus", "omicron"]
     # create dir for summary if don't exist
     if not os.path.exists(f"../../2_MD_Amber/sample/{out_dir}"):
          os.mkdir(f"../../2_MD_Amber/sample/{out_dir}")

     # define gb runs to be analyzed
     runs = ["igb2_salt150", "igb8_salt150"]

     for cmplx in complexes:
          print(f"Start PBSA analysis for {cmplx}")
          for run in runs:
               print(f"Analysing run for {run} model")
               
               # get data as a nested dictionary
               print('Checking compressed data existence....')
               mmpbsa_dir = os.path.join(input_dir, f"{cmplx}+mp1/7_mmpbsa/{run}")
               if os.path.exists(os.path.join(mmpbsa_dir, f"MMPBSA_data_{cmplx}+mp1_{run}.pickle")):
                    print("Data is compressed. Loading already pickled data...")
                    if answ == 'y':
                         pickle_mmgbsa_data(path_to_info_file=f"{mmpbsa_dir}/_MMPBSA_info",path_to_output_dir=mmpbsa_dir, cmplx=cmplx, run=run)
                    #else:data = load_pickle_data(pcike_file=f"MMPBSA_data_{cmplx}_{run}.pickle")
                    #data = load_pickle_data(pcike_file=f"{mmpbsa_dir}/MMPBSA_data_{cmplx}_{run}.pickle")
               
               else:
                    print("Data is not compressed. Compressing...")
                    pickle_mmgbsa_data(path_to_info_file=f"{mmpbsa_dir}/_MMPBSA_info",path_to_output_dir=mmpbsa_dir, cmplx=cmplx, run=run)
                    print("Compression ended")

               print("Loading data...")
               data = load_pickle_data(pcike_file=f"{mmpbsa_dir}/MMPBSA_data_{cmplx}+mp1_{run}.pickle")

               print('Finished loading data')
               
               # divide complex, receptor, ligand metrcis
               data_complex = data['gb']['complex']['TOTAL'].copy()
               data_lig = data['gb']['ligand']['TOTAL'].copy()
               data_rec = data['gb']['receptor']['TOTAL'].copy()
               # get total number of analyzed frames
               print(f"Found data for {data_complex.shape[0]} frames")

               # reshape
               data_complex = data_complex.reshape((data_complex.shape[0], 1))
               data_lig = data_lig.reshape((data_lig.shape[0], 1))
               data_rec = data_rec.reshape((data_rec.shape[0], 1))

               # construct a table
               data_table = pd.DataFrame(data=np.hstack((data_complex, data_rec, data_lig)),
                                         columns=["complex", "receptor", "ligand"])

               # calculate total binding energy
               data_table['dG'] = data_table["complex"] - data_table["receptor"] - data_table["ligand"]
               
               # assign frame number column
               data_table['frame'] = data_table.index

               # change order of columns
               data_table = data_table[['frame', 'complex', 'receptor', 'ligand', 'dG']]

               # save to csv
               data_table.to_csv(os.path.join(out_dir, f"{cmplx}+mp1_{run}.tsv"), index=False, sep='\t') # save data to tsv

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Summarize GBSA data into tsv format')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing MM-GBSA info file', required=True, default="../../2_MD_Amber/sample")

    parser.add_argument('-y', nargs=1, help='Overwrite or not already existing pickles: y/n', required=True, default="y")
    
    parser.add_argument('-o', nargs=1, help='Path to output folder', required=True, default="../tables/")
    
    args = parser.parse_args()

    sum_gbsa(input_dir=args.i[0], out_dir=args.o[0], answ=args.y[0])
