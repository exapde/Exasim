#!/bin/bash
#SBATCH --job-name=ionic2
#SBATCH --output=ionic2.out
#SBATCH --error=ionic2.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

cd app
mpirun mpiapp 1 ../datain/ ../dataout/out