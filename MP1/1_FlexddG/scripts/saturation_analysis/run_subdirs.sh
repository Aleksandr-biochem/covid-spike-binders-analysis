#!/bin/bash

for D in ./*; do
	if [ -d "$D" ]; then
		echo "Starting analysis for $D"
		cd "$D"
		sbatch flexddg.sh
		echo "Analysis for $D is finished"
		cd ..
	fi
done
