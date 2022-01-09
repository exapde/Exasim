# External packages
using Revise, DelimitedFiles, SymPy

# Add Exasim to Julia search path
cdir = pwd(); ii = findlast("Exasim", cdir);
include(cdir[1:ii[end]] * "/Installation/setpath.jl");

# Exasim packages
using Preprocessing, Mesh, Gencode, Postprocessing

# create pde structure and mesh structure
pde, mesh = Preprocessing.initializeexasim();

# Define PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";            # ModelC, ModelD, ModelW
include("pdemodel.jl");          # include the PDE model file

# Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;               # polynomial degree
pde.torder = 2;               # time-stepping order of accuracy
pde.nstage = 2;               # time-stepping number of stages
pde.dt = 0.002*ones(200);   # time step sizes
pde.soltime = collect(10:10:length(pde.dt)); # steps at which solution are collected
pde.visdt = pde.dt[1]; # visualization timestep size

gam = 1.4;                      # specific heat ratio
Re = 200;                     # Reynolds number
Pr = 0.71;                      # Prandtl number
Minf = 0.4;                     # Mach number
alpha = 5*pi/180;               # angle of attack
rinf = 1.0;                     # freestream density
ruinf = cos(alpha);             # freestream horizontal velocity
rvinf = sin(alpha);             # freestream vertical velocity
pinf = 1/(gam*Minf^2);          # freestream pressure
rEinf = 0.5+pinf/(gam-1);       # freestream energy
pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = [2.0];          # DG stabilization parameter
pde.GMRESrestart=30;            # number of GMRES restarts
pde.linearsolvertol=0.0001;     # GMRES tolerance
pde.linearsolveriter=31;        # number of GMRES iterations
pde.precMatrixType=2;           # preconditioning type
pde.NLtol = 1e-7;               # Newton tolerance
pde.NLiter=3;                   # Newton iterations

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 4;                # number of MPI processors

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

