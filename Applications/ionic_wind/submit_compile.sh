#!/bin/bash
#SBATCH --job-name=ionic_compile
#SBATCH --output=ionic_compile.out
#SBATCH --error=ionic_compile.err

source activate exasim_en
matlab -nodisplay -nodesktop -r "cd ~/Exasim/Applications/ionic_wind/; app_fullycoupled_electrons_only; save('dmd.mat','dmd'); exit"

# Usage: start in main case run directory
# sbatch submit_compile.sh
# tail -f ionic_compile.out
