"""
The script runs FlexddG output analysis in all directories
"""
import os
from subprocess import PIPE, run

# get a list of diresctories to analyze
output_dirs = [d for d in next(os.walk(f'./'))[1] if 'template' not in d]

# run analysis script in every directory
for dir in output_dirs:
	os.chdir(f'./{dir}')
	command = ['python', 'analyze_flex_ddG.py', 'output_saturation']
	result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print(f'Directory {dir} analyzed')
	os.chdir('..')
