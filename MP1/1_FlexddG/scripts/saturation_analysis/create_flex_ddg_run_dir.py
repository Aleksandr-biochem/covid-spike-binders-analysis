from pyxmolpp2 import PdbFile, rId
from subprocess import call
import os
import fileinput
import sys

# Generate function to iterate over multiple inputs (simultaniously create subdirs for several inputs)
## each iteration - print which experiment is now created
def iterate_inputs(file, search_exp, replace_exp):
    for line in fileinput.input(file, inplace=1):
        if search_exp in line:
            line = line.replace(search_exp, replace_exp)
        # write for which experiment now directories are created
        sys.stdout.write(line)

# Function to create subdirs for all residues
def create_dir(parent, daughter):
    call(['cp', '-R', parent, daughter])

path_to_pdb = "template/inputs/pdb_deltaplus/delta_plus+mp1_minimized.pdb"
pdb = PdbFile(path_to_pdb).frames()[0]

path_to_template = "template/."
# get only first 36 residues from A chain - they're mutated
interest_res = pdb.residues.filter(rId.is_in(set(range(1, 37))))
for res in interest_res:
    out_dir = f'{res.id.serial:02d}_{res.name}'
    # recursively create subfolder of  curr residue in curr dir
    os.makedirs(out_dir, exist_ok=True)
        
    create_dir(path_to_template, out_dir)
        
    # make copies of scripts
    run_script = os.path.join(out_dir, "run_saturation.py")
    example_mutation_string = "residue_to_mutate = ('A', 29, '')"
    current_mutation_string = f"residue_to_mutate = ('A', {res.id.serial}, '')"
    iterate_inputs(run_script, example_mutation_string, current_mutation_string)

    bash_script = os.path.join(out_dir, "flexddg.sh")
    example_job_name = '#SBATCH --job-name="flex_ddg_S29"'
    current_job_name = f'#SBATCH --job-name="{res.id.serial:02d}_{res.name}"'
    iterate_inputs(bash_script, example_job_name, current_job_name)

    plot_script = os.path.join(out_dir, "plot_dG.py")
    example_res_id = "residue_number = 29"
    current_res_id = f"residue_number = {res.id.serial:02d}"
    iterate_inputs(plot_script, example_res_id, current_res_id)

    example_res_name = 'initial_residue = "SER"'
    current_res_name = f'initial_residue = "{res.name}"'
    iterate_inputs(plot_script, example_res_name, current_res_name) 
