#!/bin/bash
# Using the array job feature of SLURM to submit a serial farm.
# The input file table.dat contains individual cases - one case per line
# A simple test table can be generated with
# for ((i=1; i<=4; i++)); do echo "sleep $(($RANDOM % 30))"; done > table.dat

export TABLE=./table.dat

# Total number of cases (= number of jobs to submit):
N_cases=$(cat "$TABLE" | wc -l)

# Submitting an array job to the scheduler:
sbatch --array=1-$N_cases ./run_many.sh