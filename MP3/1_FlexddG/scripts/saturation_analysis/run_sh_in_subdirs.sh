#!/bin/bash

for D in ./*; do
    if [ -d "$D" ]; then
        cd "$D"
        sbatch flexddg.sh
        cd ..
    fi
done
