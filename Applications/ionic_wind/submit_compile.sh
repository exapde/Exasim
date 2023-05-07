#!/bin/bash
#SBATCH --job-name=ionic
#SBATCH --output=ionic.out
#SBATCH --error=ionic.err

source activate exasim_en
matlab -nodisplay -nodesktop -r "cd ~/Exasim/Applications/ionic_wind_fully_coupled/; pde_fullycoupled_3eqns; save('dmd.mat','dmd'); exit"