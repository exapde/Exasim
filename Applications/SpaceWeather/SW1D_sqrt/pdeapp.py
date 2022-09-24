"""
Module that runs GITM 1D (reduced stiffness formulation) model.

Authors: Opal Issan
Version: Sept 9, 2022
"""
# import external modules
import os, sympy
from numpy import *
from pdeparams import pdeparams

# Add Exasim to Python search path
#cdir = '/Users/oissan/PycharmProjects/Exasim/Applications/SpaceWeather/SW1D_sqrt';
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii + 6)] + "/Installation/setpath.py").read())

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde, mesh = Preprocessing.initializeexasim()

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelD"  # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"  # name of a file defining the PDE model

# Choose computing platform and set number of processors
pde['platform'] = "cpu"
pde['mpiprocs'] = 1  # number of MPI processors

pde, mesh = pdeparams(pde,mesh)

# search compilers and set options
pde = Gencode.setcompilers(pde);

# generate input files and store them in datain folder
pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

# generate source codes and store them in app folder
Gencode.gencode(pde);

# compile source codes to build an executable file and store it in app folder
compilerstr = Gencode.compilecode(pde);

# run executable file to compute solution and store it in dataout folder
runstr = Gencode.runcode(pde,1);

# get solution from output files in dataout folder
sol = Postprocessing.fetchsolution(pde,master,dmd,"dataout");
