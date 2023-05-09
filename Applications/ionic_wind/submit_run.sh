#!/bin/bash
#SBATCH --job-name=ionic_run
#SBATCH --output=ionic_run.out
#SBATCH --error=ionic_run.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10

cd app
mpirun mpiapp 10 ./datain/ ./dataout/out

# Usage: start in main case run directory
# sbatch submit_run.sh
# tail -f ionic_run.out