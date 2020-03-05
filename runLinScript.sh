#!/bin/bash

#SBATCH --time=05:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH --qos=pws

python3 runLinearRegressions.py $1 $2