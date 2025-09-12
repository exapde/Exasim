# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/install/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim();

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelC";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 2;          # polynomial degree
pde.torder = 1;          # time-stepping order of accuracy
pde.nstage = 1;          # time-stepping number of stages
pde.dt = 0.003*ones(10000);   # time step sizes
pde.saveSolFreq = 10;          # solution is saved every 10 time steps
pde.soltime = collect(10:10:length(pde.dt)); # steps at which solution are collected
pde.visdt = 0.003; # visualization timestep size

gam = 1.4;              # specific heat ratio
Minf = 0.5;             # Mach number
rinf = 1.0;             # freestream density
uinf = 1.0;             # freestream horizontal velocity
vinf = 0.0;             # freestream vertical velocity
pinf = 1/(gam*Minf^2);  # freestream pressure
rEinf = 0.5+pinf/(gam-1); # freestream energy
pde.physicsparam = [gam Minf rinf uinf vinf rEinf];
pde.tau = [2.0];          # DG stabilization parameter
pde.GMRESrestart=30;            # number of GMRES restarts
pde.linearsolvertol=0.0001;     # GMRES tolerance
pde.linearsolveriter=31;        # number of GMRES iterations
pde.precMatrixType=2;           # preconditioning type
pde.NLtol = 1e-7;               # Newton tolerance
pde.NLiter=3;                   # Newton iterations

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# read a grid from a file
mesh.p,mesh.t = Mesh.readmesh("grid.bin",0);
# expressions for domain boundaries
mesh.boundaryexpr = [p -> (sqrt.((p[1,:].-0.5).*(p[1,:].-0.5)+(p[2,:].*p[2,:])) .< 3), p -> (p[1,:] .< 20)];
mesh.boundarycondition = [1 2];
# expressions for curved boundaries
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = [p -> p[2,:].^2-(5*0.01*12*(0.29690*sqrt.(abs.(p[1,:]))-0.12600*p[1,:]-0.35160*p[1,:].^2+0.28430*p[1,:].^3-0.10150*p[1,:].^4)).^2, p -> 0];

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
pde.visscalars = ["density", 1, "energy", 4];  # list of scalar fields for visualization
pde.visvectors = ["momentum", [2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
