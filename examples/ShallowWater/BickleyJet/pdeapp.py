# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/install/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim();

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelC";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 1;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 4;         # polynomial degree
pde['torder'] = 2;          # time-stepping order of accuracy
pde['nstage'] = 2;          # time-stepping number of stages
pde['dt'] = 0.02*numpy.ones(5000);   # time step sizes
pde['soltime'] = numpy.arange(50,pde['dt'].size+1,50); # steps at which solution are collected
pde['visdt'] = 1.0; # visualization timestep size

gam = 10.0;                      # gravity
pde['physicsparam'] = [gam];
pde['tau'] = numpy.array([1.0]); # DG stabilization parameter
pde['GMRESrestart']=15;            # number of GMRES restarts
pde['linearsolvertol']=1e-12;      # GMRES tolerance
pde['linearsolveriter']=16;        # number of GMRES iterations
pde['precMatrixType']=2;           # preconditioning type
pde['NLtol'] = 1e-12;              # Newton tolerance
pde['NLiter']=2;                   # Newton iterations

# create a mesh of 10 by 10 quads on a square domain
mesh['p'], mesh['t'] = Mesh.SquareMesh(64,64,1)[0:2];
pi = numpy.pi;
mesh['p'] = (4*pi)*mesh['p'] - 2*pi;
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < -2*pi+1e-3), lambda p: (p[0,:] > 2*pi-1e-3), lambda p: (p[1,:] > 2*pi-1e-3), lambda p: (p[0,:] < -2*pi+1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary
mesh['periodicexpr'] = [[2, lambda p: p[1,:], 4, lambda p: p[1,:]], [1, lambda p: p[0,:], 3, lambda p: p[0,:]]];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["density", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["momentum", [1, 2]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");

