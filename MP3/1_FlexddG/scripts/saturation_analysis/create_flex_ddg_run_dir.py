"""
The script creates subdirectories for saturation analysis with FlexddG using template folder with scripts
"""
from pyxmolpp2 import PdbFile, rId
from subprocess import call
import os
import fileinput
import sys


def replace_line(file, search_exp, replace_exp):
    for line in fileinput.input(file, inplace=1):
        if search_exp in line:
            line = line.replace(search_exp, replace_exp)
        sys.stdout.write(line)


def cp_dir(source, target):
    call(['cp', '-R', source, target])

path_to_pdb = "template/inputs/pdb_deltaplus/delta_plus+mp3_minimized.pdb" # path to structure for saturation analysis
path_to_template = "template/." # path to template files

pdb = PdbFile(path_to_pdb).frames()[0] # load pdb structure

# here ids for residues of MP3 helices on interface with RBD are set manually as 1-41
res_of_interest = pdb.residues.filter(rId.is_in(set(range(1, 42))))

for res in res_of_interest:
    # a dir is created for each residue mutation 
    out_dir = f"{res.id.serial:02d}_{res.name}"
    os.makedirs(out_dir, exist_ok=True)
    cp_dir(path_to_template, out_dir) # copy all files to residue folder

    run_file = os.path.join(out_dir, "run_saturation.py")
    example_mutation_string = "residue_to_mutate = ('A', 29, '')"
    current_mutation_string = f"residue_to_mutate = ('A', {res.id.serial}, '')"
    # the example line in run_saturation.py specifying mutation position is modifyed according to current residue
    replace_line(run_file, example_mutation_string, current_mutation_string)
    
    sbatch_file = os.path.join(out_dir, "flexddg.sh")
    example_job_name = '#SBATCH --job-name="flex_ddg_S29"'
    current_job_name = f'#SBATCH --job-name="{res.id.serial:02d}_{res.name}"'
    # the example line in flexddg.sh specifying job name is modifyed according to current residue
    replace_line(sbatch_file, example_job_name, current_job_name)
    
    plot_file = os.path.join(out_dir, "plot_dG.py")
    example_res_id = "residue_number = 29"
    current_res_id = f"residue_number = {res.id.serial:02d}"
    # the example line with res id in plot_dG.py is modifyed according to current residue
    replace_line(plot_file, example_res_id, current_res_id)

    example_res_name = 'initial_residue = "SER"'
    current_res_name = f'initial_residue = "{res.name}"'
    # the example line with res name in plot_dG.py is modifyed according to current residue
    replace_line(plot_file, example_res_name, current_res_name)

