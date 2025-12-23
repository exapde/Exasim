#!/usr/bin/env bash

set -e  # stop on first error
set -u  # error on unset variables

echo "Running Poisson 1D..."
../build/text2code poisson/poisson1d/pdeapp.txt
../build/cput2cEXASIM poisson/poisson1d/pdeapp.txt | tee poisson1d.log

echo "Running Poisson 2D..."
../build/text2code poisson/poisson2d/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM poisson/poisson2d/pdeapp.txt | tee poisson2d.log

echo "Running Poisson 3D..."
../build/text2code poisson/poisson3d/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM poisson/poisson3d/pdeapp.txt | tee poisson3d.log

echo "Running Poisson Periodic..."
../build/text2code poisson/periodic/pdeapp.txt
mpirun -np 2 ../build/cpumpit2cEXASIM poisson/periodic/pdeapp.txt | tee poissonperiodic.log

echo "Running Poisson orion..."
../build/text2code poisson/orion/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM poisson/orion/pdeapp.txt | tee poissonorion.log

echo "Running Poisson lshape..."
../build/text2code poisson/lshape/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM poisson/lshape/pdeapp.txt | tee poissonlshape.log

echo "Running Poisson cone..."
../build/text2code poisson/cone/pdeapp.txt
mpirun -np 8 ../build/cpumpit2cEXASIM poisson/cone/pdeapp.txt | tee poissoncone.log

echo "Running Navier-Stokes orion..."
../build/text2code navierstokes/orion/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM navierstokes/orion/pdeapp.txt | tee nsorion.log

echo "Running Navier-Stokes nsmach8..."
../build/text2code navierstokes/nsmach8/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM navierstokes/nsmach8/pdeapp.txt | tee nsmach8.log

echo "Running Navier-Stokes naca0012steady..."
../build/text2code navierstokes/naca0012steady/pdeapp.txt
mpirun -np 4 ../build/cpumpit2cEXASIM navierstokes/naca0012steady/pdeapp.txt | tee naca0012steady.log

echo "Running Navier-Stokes naca0012unsteady..."
../build/text2code navierstokes/naca0012unsteady/pdeapp.txt
mpirun -np 8 ../build/cpumpit2cEXASIM navierstokes/naca0012unsteady/pdeapp.txt | tee naca0012unsteady.log

echo "All runs completed successfully."