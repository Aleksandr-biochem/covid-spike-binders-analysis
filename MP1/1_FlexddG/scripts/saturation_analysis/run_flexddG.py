import os
from subprocess import PIPE, run

# get a list of diresctories to analyze
output_dirs = [d for d in next(os.walk(f'./'))[1] if 'template' not in d]

# run analysis script in every directory
for directory in output_dirs:
	os.chdir(f'./{directory}')
	command = ['python3', 'analyze_flex_ddG.py', 'output_saturation']
	print(f'Saturation analysis for {directory} has started')
	result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	# state of analysis
	print(f'Directory {directory} is analyzed')
	os.chdir('..')