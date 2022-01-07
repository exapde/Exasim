# Specify an Exasim version to run
version = "Version0.1";

# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/Installation/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create pde object and mesh object
pde,mesh = Preprocessing.initializeexasim(version);

# Define a PDE model: governing equations and boundary conditions
pde['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 4;        # number of MPI processors

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['torder'] = 2;          # time-stepping order of accuracy
pde['nstage'] = 2;          # time-stepping number of stages
pde['dt'] = 0.002*numpy.ones(10000);   # time step sizes
pde['soltime'] = numpy.arange(10,pde['dt'].size+1,10); # steps at which solution are collected
pde['visdt'] = pde['dt'][0]; # visualization timestep size

gam = 1.4;                      # specific heat ratio
Re = 10000;                     # Reynolds number
Pr = 0.71;                      # Prandtl number
Minf = 0.4;                     # Mach number
alpha = 5*numpy.pi/180;               # angle of attack
rinf = 1.0;                     # freestream density
ruinf = numpy.cos(alpha);             # freestream horizontal velocity
rvinf = numpy.sin(alpha);             # freestream vertical velocity
pinf = 1/(gam*Minf**2);          # freestream pressure
rEinf = 0.5+pinf/(gam-1);       # freestream energy
pde['physicsparam'] = [gam, Re, Pr, Minf, rinf, ruinf, rvinf, rEinf];
pde['tau'] = numpy.array([2.0]); # DG stabilization parameter
pde['GMRESrestart']=30;            # number of GMRES restarts
pde['linearsolvertol']=0.0001;     # GMRES tolerance
pde['linearsolveriter']=31;        # number of GMRES iterations
pde['precMatrixType']=2;           # preconditioning type
pde['NLtol'] = 1e-7;               # Newton tolerance
pde['NLiter']=3;                   # Newton iterations

# read a grid from a file
mesh['p'], mesh['t'] = Mesh.readmesh("grid.bin",0);
mesh['t'] = mesh['t'] - 1;

# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (numpy.sqrt((p[0,:]-0.5)*(p[0,:]-0.5)+p[1,:]*p[1,:]) < 3), lambda p: (p[0,:] < 20)];
mesh['boundarycondition'] = numpy.array([1, 2]); # Set boundary condition for each boundary
#expressions for curved boundaries
mesh['curvedboundary'] = numpy.array([1, 0]);
mesh['curvedboundaryexpr'] = [lambda p: p[1,:]**2-(5*0.01*12*(0.29690*numpy.sqrt(numpy.abs(p[0,:]))-0.12600*p[0,:]-0.35160*p[0,:]**2+0.28430*p[0,:]**3-0.10150*p[0,:]**4))**2, lambda p: 0];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["density", 0, "energy", 3]; # list of scalar fields for visualization
pde['visvectors'] = ["momentum", [1, 2]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
