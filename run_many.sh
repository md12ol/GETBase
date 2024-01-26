#!/bin/bash
# Here you should provide the sbatch arguments to be used in all jobs in this serial farm
#SBATCH --time=30:00:00
#SBATCH --mem=1G
#SBATCH --account=def-houghten
#SBATCH --cpus-per-task=1

# Case number is given by $SLURM_ARRAY_TASK_ID environment variable:
i=$SLURM_ARRAY_TASK_ID

# Extracing the $i-th line from file $TABLE:
LINE=`sed -n ${i}p "$TABLE"`
# Echoing the command (optional), with the case number prepended:
echo "$i; $LINE"
# Executing the command:
eval "$LINE"