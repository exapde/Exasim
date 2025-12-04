#!/bin/bash

# export CASE_PREFIX=009
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX')"

# Case just to double check that the source term is still the way it was
# export CASE_PREFIX=010
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX')"

# First attempt at sponge layer
# export CASE_PREFIX=011
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Flipped the sign on beta
# export CASE_PREFIX=012
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Beta too large?
# export CASE_PREFIX=013
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Crash on symmetry boundary, going to visualize the sponge layer
# export CASE_PREFIX=014
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 20)"

# # Crash on symmetry boundary, going to visualize all sponge source terms
# export CASE_PREFIX=015
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Forget it, going to add several more solution fields to diagnose why the sponge layer is crashing
# export CASE_PREFIX=016
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 10)"

# # Forget it, going to add several more solution fields to diagnose why the sponge layer is crashing
# export CASE_PREFIX=017
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # That is working, but it looks like the source is ADDING to the wave, not subtracting. Let's try flipping the sign.
# export CASE_PREFIX=018
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Proper sign, but increase sponge strength and change to a radial filter instead of purely in x bc of wall reflections
# export CASE_PREFIX=018
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # First test of the tdep energy source
# export CASE_PREFIX=019
# matlab -batch "pdeapp_small(1, 0, 0, 1e-1, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Crash almost immediately, try reducing the timestep
# export CASE_PREFIX=020
# matlab -batch "pdeapp_small(1, 0, 0, 1e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Source term sign should be flipped -- energy was decreasing
# export CASE_PREFIX=021
# matlab -batch "pdeapp_small(1, 0, 0, 1e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # That worked, now with a much higher energy deposition rate
# export CASE_PREFIX=022
# matlab -batch "pdeapp_small(1, 0, 0, 1e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Writing more diagnostic fields to the output file
# export CASE_PREFIX=023
# matlab -batch "pdeapp_small(1, 0, 0, 1e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# # Tightening up sigmoid slope and decreasing kernel radius to match the paper
# export CASE_PREFIX=024
# matlab -batch "pdeapp_small(1, 0, 0, 1e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 1000)"

# Stat of debugging negative density
export CASE_PREFIX=025
matlab -batch "pdeapp_small(1, 0, 0, 0,5e-2, '/home/austinsp/code112425/Exasim/examples/euler_electrode/bamg_371k.msh', 1, '/home/austinsp/code112425/Exasim/build/$CASE_PREFIX', 300)"


sbatch --export=ALL 2_submit.sh
