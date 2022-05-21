import os
import argparse
import pandas as pd

def extract_data(path, dir, state, backrup, n_struct, score_func):
    
    scores, aa = dict(), 'GAVLISTMCNQDEKRHFWYP'
    for a in aa:
        scores[a] = 0

    def add_value(aa, ddg):
        scores[aa] = ddg    
    
    analysis_dir = [d for d in next(os.walk(f'{path}/{dir}'))[1] if 'analysis_output' in d]

    if len(analysis_dir) > 0:
        df = pd.read_csv(f'{path}/{dir}/{analysis_dir[0]}/output_saturation-results.csv')
        
        display_columns = ['backrub_steps', 'case_name', 'nstruct', 'score_function_name', 'scored_state', 'total_score']
        ddg_data = df.loc[ df['scored_state'] == state ][display_columns]
        
        filtered = ddg_data[(ddg_data['backrub_steps'] == backrup)
                            & (ddg_data['nstruct'] == n_struct)
                            & (ddg_data['score_function_name'] == score_func)]

        filtered.apply(lambda row: add_value(row.case_name[-1], row.total_score), axis=1)

        index = int(dir[:2])
        residue = f'{dir[-3]}{dir[-2:].lower()}'
        return [index, residue] + [float(item[1]) for item in scores.items()]
    else:
        print(f'Directory {dir} contains no analysis output')

def retrieve_mutation_data(path, state, backrup, n_struct, score_func):
    
    aa_dict = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q',
               'Lys': 'K', 'Ile': 'I', 'Pro': 'P', 'Thr': 'T',
               'Phe': 'F', 'Asn': 'N', 'Gly': 'G', 'His': 'H',
               'Leu': 'L', 'Arg': 'R', 'Trp': 'W', 'Ala': 'A',
               'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

    sub_dirs = next(os.walk(path))[1] # list subfolders with data
    data = []
    for dir in sub_dirs:
        if dir != 'template':
            new_row = extract_data(path, dir, state, backrup, n_struct, score_func)
            if new_row is not None:
                new_row[1] = aa_dict[new_row[1]] # change to single-letter code
                data.append(new_row)
    
    residues =  ['G', 'A', 'V', 'L', 'I', # aliphatic radical
                 'S', 'T',               # oh- radical
                 'M', 'C',               # containing S
                 'N', 'Q',               # carboxamide
                 'D', 'E',               # acidic 
                 'K', 'R', 'H',          # basic
                 'F', 'W', 'Y',          # aromatic
                 'P']                    # imine

    cols =  ['index', 'residue'] + residues
    data = pd.DataFrame(data=data, columns=cols)
    return data

if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(description='Filter and summarize information from FlexddG output folders:')
    
    parser.add_argument('-i', nargs=1, help='Path to folder containing subfolders with FlexddG output', required=True)
    
    parser.add_argument('-s', nargs=1, help='State to summarize: mut_dG, wt_dG or ddG', required=True)
    
    parser.add_argument('-b', nargs=1,   help='Value for backrub_steps to filter', required=True)
    
    parser.add_argument('-n', nargs=1,   help='Values for nstruct to filter', required=True)
    
    parser.add_argument('-f', nargs=1,   help='Value for score_function to filter: fa_talaris2014 or fa_talaris2014-gam', required=True)
    
    parser.add_argument('-o', nargs=1,   help='Output .csv name', required=True)
    
    args = parser.parse_args()
    mutation_data = retrieve_mutation_data(args.i[0], args.s[0], int(args.b[0]), int(args.n[0]), args.f[0])
    mutation_data.to_csv(f'{args.o[0]}.csv')
