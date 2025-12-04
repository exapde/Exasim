#!/bin/bash
#SBATCH --job-name=exasim
#SBATCH --output=logs/run_%A_%a.out  # %A=JobArrayID, %a=Index
#SBATCH --error=logs/run_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:h200:1

module load openblas/0.3.26

# --- CONFIGURATION ---
BUILD_DIR="/home/austinsp/code112425/Exasim/build/$CASE_PREFIX"

cd "$BUILD_DIR"

# Build Step
cmake -D EXASIM_NOMPI=ON \
        -D EXASIM_MPI=OFF \
        -D EXASIM_CUDA=ON \
        ../../install && \
cmake --build .

./gpuEXASIM 1 "$BUILD_DIR/datain/" "$BUILD_DIR/dataout/out"
