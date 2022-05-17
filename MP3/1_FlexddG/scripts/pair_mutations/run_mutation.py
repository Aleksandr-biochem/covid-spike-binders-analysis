#!/usr/bin/python

from __future__ import print_function

import socket
import sys
import os
import subprocess

use_multiprocessing = True
if use_multiprocessing:
    import multiprocessing
    max_cpus = 20 # We might want to not run on the full number of cores, as Rosetta take about 2 Gb of memory per instance

###################################################################################################################################################################
# Important: The variables below are set to values that will make the run complete faster (as a tutorial example), but will not give scientifically valid results.
#            Please change them to the "normal" default values before a real run.
###################################################################################################################################################################

rosetta_scripts_path = os.path.expanduser("/home/olebedenko/tools/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/rosetta_scripts.mpi.linuxgccrelease")
nstruct = 35 # Normally 35
max_minimization_iter = 5000 # Normally 5000
abs_score_convergence_thresh = 1.0 # Normally 1.0
number_backrub_trials = 35000 # Normally 35000
backrub_trajectory_stride = 17500 # Can be whatever you want, if you would like to see results from earlier time points in the backrub trajectory. 7000 is a reasonable number, to give you three checkpoints for a 35000 step run, but you could also set it to 35000 for quickest run time (as the final minimization and packing steps will only need to be run one time).
path_to_script = 'ddG-backrub.xml'


mutation_positions = {3: 'D', 7: 'M', 10: 'T', 17: 'L', 34: 'E', 37: 'D'} # residues to be mutated

# Residue combinations to mutate
# Format specification: a list of (Chain, tuple(PDB residue numbers), tuple(new residues in single letter notation).
residue_combinations_to_mutate = [('A', (37, 3), ('R', 'W')), ('A', (37, 3), ('R', 'H')),
                                  ('A', (37, 7), ('R', 'M')), ('A', (37, 7), ('R', 'F')),
                                  ('A', (37, 10), ('R', 'Y')), ('A', (37, 10), ('R', 'W')), ('A', (37, 10), ('R', 'F')),
                                  ('A', (37, 10), ('R', 'H')), ('A', (37, 10), ('R', 'R')), ('A', (37, 10), ('R', 'Q')),
                                  ('A', (37, 17), ('R', 'R')), ('A', (37, 34), ('R', 'W'))]

if not os.path.isfile(rosetta_scripts_path):
    print('ERROR: "rosetta_scripts_path" variable must be set to the location of the "rosetta_scripts" binary executable')
    print('This file might look something like: "rosetta_scripts.linuxgccrelease"')
    raise Exception('Rosetta scripts missing')

def run_flex_ddg(name, input_path, input_pdb_path, chains_to_move, mut_aa, nstruct_i ):

    mutation_residues, mutation_icodes =  mut_aa # unpack residues to insert
    res1, res2 = mutation_residues # unpack residues
    code1, code2 = mutation_icodes # unpack positions

    output_directory = os.path.join( 'output_saturation', os.path.join( '%s_%d%s' % (name, code2, res2), '%02d' % nstruct_i ) )
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    # positions of mutation in initial protein
    mutation_chain = 'A'
    mutation_res1, mutation_res2 = mutation_positions[code1], mutation_positions[code2]

    resfile_path = os.path.join( output_directory, 'mutate_%s%d%s_to_%s.resfile' % (mutation_chain, code2, mutation_res2, res2) )
    with open( resfile_path, 'w') as f:
        f.write( 'NATRO\nstart\n')
        f.write( '%d %s PIKAA %s\n' % (code1, mutation_chain, res1) )
        f.write( '%d %s PIKAA %s\n' % (code2, mutation_chain, res2) )

    flex_ddg_args = [
        os.path.abspath(rosetta_scripts_path),
        "-s %s" % os.path.abspath(input_pdb_path),
        '-parser:protocol', os.path.abspath(path_to_script),
        '-parser:script_vars',
        'chainstomove=' + chains_to_move,
        'mutate_resfile_relpath=' + os.path.abspath( resfile_path ),
        'number_backrub_trials=%d' % number_backrub_trials,
        'max_minimization_iter=%d' % max_minimization_iter,
        'abs_score_convergence_thresh=%.1f' % abs_score_convergence_thresh,
        'backrub_trajectory_stride=%d' % backrub_trajectory_stride ,
        '-restore_talaris_behavior',
        '-in:file:fullatom',
        '-ignore_unrecognized_res',
        '-ignore_zero_occupancy false',
        '-ex1',
        '-ex2',
    ]

    log_path = os.path.join(output_directory, 'rosetta.out')

    print( 'Running Rosetta with args:' )
    print( ' '.join(flex_ddg_args) )
    print( 'Output logged to:', os.path.abspath(log_path) )
    print()

    outfile = open(log_path, 'w')
    process = subprocess.Popen(flex_ddg_args, stdout=outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = output_directory)
    returncode = process.wait()
    outfile.close()

if __name__ == '__main__':
    
    cases = []
    for nstruct_i in range(1, nstruct + 1 ):
        for case_name in os.listdir('inputs'):
            case_path = os.path.join( 'inputs', case_name )
            for f in os.listdir(case_path):
                if f.endswith('.pdb'):
                    input_pdb_path = os.path.join( case_path, f )
                    break

            with open( os.path.join( case_path, 'chains_to_move.txt' ), 'r' ) as f:
                chains_to_move = f.readlines()[0].strip()

            for combination in residue_combinations_to_mutate:
                mutation_chain, mutation_icodes, mutation_residues = combination 
                res1, res2 = mutation_residues # unpack pairs residues
                code1, code2 = mutation_icodes # unpack pairs positions
                cases.append( ('%s_%s%d%s' % (case_name, mutation_chain, code1, res1), case_path, input_pdb_path, chains_to_move, (mutation_residues, mutation_icodes), nstruct_i) )

    if use_multiprocessing:
        pool = multiprocessing.Pool( processes = min(max_cpus, multiprocessing.cpu_count()) )

    for args in cases:
        if use_multiprocessing:
            pool.apply_async( run_flex_ddg, args = args )
        else:
            run_flex_ddg( *args )

    if use_multiprocessing:
        pool.close()
        pool.join()
