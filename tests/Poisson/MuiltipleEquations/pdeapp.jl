# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/Installation/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

pde = Array{Any, 1}(undef, 2);
mesh = Array{Any, 1}(undef, 2);

# create pde and mesh for each PDE model
include("pdeapp1.jl"); 
include("pdeapp2.jl"); 

# call exasim to generate and run C++ code to solve the PDE models
sol,pde,mesh,master,dmd,compilerstr,runstr = exasim(pde,mesh);

