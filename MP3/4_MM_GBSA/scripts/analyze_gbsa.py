from MMPBSA_mods import API as MMPBSA_API
import pandas as pd
import numpy as np
import pickle
import os

# a function to convert data in _MMPBSA_info to pickle for more rapid parsing
def pickle_mmgbsa_data(path_to_info_file, path_to_output_dir=".", outname="MMPBSA_data.pickle"):
    wd = os.getcwd()
    os.chdir(os.path.dirname(path_to_info_file))
    data = MMPBSA_API.load_mmpbsa_info("_MMPBSA_info")
    os.chdir(wd)
    os.makedirs(path_to_output_dir, exist_ok=True)
    with open(os.path.join(path_to_output_dir, outname), 'wb') as f:
        pickle.dump(data, f)

# function to load pickle data
def load_pickle_data(pcike_file="MMPBSA_data.pickle"):
    with open(pcike_file, 'rb') as f:
        data = pickle.load(f)
    return data

# define folder-mutants names to be analyzed 
mutants = ["wt+mp3", "alpha+mp3", "delta+mp3", "delta_plus+mp3", "omicron+mp3", "delta_p+mp3_d37r", "delta_p+mp3_d37r_t10w"]

# define gb runs to be analyzed
runs = ["igb2_salt150", "igb8_salt150"]

os.makedirs('../mmgbsa_output', exist_ok=True)

# analyze each mutant and each run
for mutant in mutants:
     print(f"Analysing {mutant}")
     for run in runs:
          print(f"Analysing run {run}")

          # change working dir to the analysis folder
          # os.chdir(f"../../2_MD_Amber/{mutant}/7_mmpbsa/{run}")

          # dump data into pickle if no pickle is found
          if not os.path.isfile(f"../../2_MD_Amber/{mutant}/7_mmpbsa/{run}/MMPBSA_data_{mutant}_{run}.pickle"):
               print("Dumping data to pickle")
               pickle_mmgbsa_data(path_to_info_file=f"../../2_MD_Amber/{mutant}/7_mmpbsa/{run}/_MMPBSA_info",
                                  path_to_output_dir=f"../../2_MD_Amber/{mutant}/7_mmpbsa/{run}",
                                  outname=f"MMPBSA_data_{mutant}_{run}.pickle")

          # load decomp_data
          data = load_pickle_data(pcike_file=f"../../2_MD_Amber/{mutant}/7_mmpbsa/{run}/MMPBSA_data_{mutant}_{run}.pickle")

          data_complex = data['gb']['complex']['TOTAL'].copy() # total contributions data for complex
          data_lig = data['gb']['ligand']['TOTAL'].copy() # total contributions data for ligand
          data_rec = data['gb']['receptor']['TOTAL'].copy() # total contributions data for receptor
          print(f"Found data for {data_complex.shape[0]} frames")

          # reshape data
          data_complex = data_complex.reshape((data_complex.shape[0], 1))
          data_lig = data_lig.reshape((data_lig.shape[0], 1))
          data_rec = data_rec.reshape((data_rec.shape[0], 1))

          # construct a table of 
          data_table = pd.DataFrame(data=np.hstack((data_complex, data_rec, data_lig)),
                                    columns=["complex", "receptor", "ligand"])

          # calculate total binding energy
          data_table['delta'] = data_table["complex"] - (data_table["receptor"] + data_table["ligand"])

          # save data to .csv in initial dir
          data_table.to_csv(f"../mmgbsa_output/{mutant}_{run}.csv")

