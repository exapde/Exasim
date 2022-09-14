"""
Module that runs GITM 1D (reduced stiffness formulation) model.

Authors: Opal Issan
Version: Sept 9, 2022
"""
# import external modules
import numpy as np
import os, sympy
from mesh1D_adapted import mesh1D_adapted

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

# Set discretization parameters, physical parameters, and solver parameters
pde["porder"] = 2  # polynomial degree
pde["torder"] = 2  # time - stepping order of  accuracy
pde["nstage"] = 2  # time - stepping number of stages
pde["dt"] = 0.190116193978700 * np.ones(86390)  # time step sizes
pde["visdt"] = 0.190116193978700  # visualization timestep size
pde["saveSolFreq"] = 360  # solution is saved every 100 time steps
pde["soltime"] = 360 * np.arange(1, 241)  # steps at which solution are collected
pde["timestepOffset"] = 0

# Vectors of physical and external parameters
# [gam Gr Pr Fr Keuv  M  rho T0 T1 R0 R1 H EUVeff model longitude latitude declinationSun tauA t0];
pde["physicsparam"] = np.loadtxt('inputs/physparam.csv', delimiter=',')

# External params (EUV)  [lambda,crossSections,AFAC,F74113]
pde["externalparam"] = np.loadtxt('inputs/externalparam.csv', delimiter=',')

print(pde["externalparam"])

# Solver parameters
pde["extStab"] = 1
pde["tau"] = 1.0  # DG stabilization parameter
pde["GMRESrestart"] = 29  # number of GMRES restarts
pde["linearsolvertol"] = 1e-16  # GMRES tolerance
pde["linearsolveriter"] = 30  # number of GMRES iterations
pde["precMatrixType"] = 2  # preconditioning type
pde["NLtol"] = 1e-10  # Newton toleranccd dataoue
pde["NLiter"] = 2  # Newton iterations
pde["matvectol"] = 1e-7
pde["RBdim"] = 8

R0, R1 = 5.918775593756355 * 1e2, 6.375605188371779 * 1e2
# need to debug this.
mesh['p'],mesh['t'] = mesh1D_adapted(R0, R1)
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[0,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3)]
mesh['boundarycondition'] = np.array([1, 2])  # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh = Postprocessing.exasim(pde, mesh)[0:3]
