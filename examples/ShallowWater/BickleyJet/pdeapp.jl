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
pde.porder = 4;          # polynomial degree
pde.torder = 2;          # time-stepping order of accuracy
pde.nstage = 2;          # time-stepping number of stages
pde.dt = 0.02*ones(5000);   # time step sizes
pde.soltime = collect(50:50:length(pde.dt)); # steps at which solution are collected
pde.visdt = 1; # visualization timestep size

gam = 10.0;                      # gravity
pde.physicsparam = [gam 0.0];
pde.tau = [1.0];          # DG stabilization parameter
pde.GMRESrestart=15;            # number of GMRES restarts
pde.linearsolvertol=1e-12;      # GMRES tolerance
pde.linearsolveriter=16;        # number of GMRES iterations
pde.precMatrixType=2;           # preconditioning type
pde.NLtol = 1e-12;              # Newton tolerance
pde.NLiter=2;                   # Newton iterations

# Choose computing platform and set number of processors
#pde.platform = "gpu";           # choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;                # number of MPI processors

# create a linear mesh for a square domain
mesh.p, mesh.t = Mesh.SquareMesh(64,64,1); # a mesh of 8 by 8 quadrilaterals
mesh.p = (4*pi)*mesh.p .- 2*pi;
# expressions for disjoint boundaries
mesh.boundaryexpr = [p -> (p[2,:] .< -2*pi+1e-3), p -> (p[1,:] .> 2*pi-1e-3), p -> (p[2,:] .> 2*pi-1e-3), p -> (p[1,:] .< -2*pi+1e-3)];
mesh.boundarycondition = [1 1 1 1]; # Set boundary condition for each disjoint boundary
mesh.periodicexpr = [2 p->p[2,:] 4 p->p[2,:]; 1 p->p[1,:] 3 p->p[1,:]];

# call exasim to generate and run C++ code to solve the PDE model
#sol, pde, mesh,~,~,~,~  = Postprocessing.exasim(pde,mesh);

# search compilers and set options
pde = Gencode.setcompilers(pde);

# generate input files and store them in datain folder
pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh);

# generate source codes and store them in app folder
Gencode.gencode(pde);

# compile source codes to build an executable file and store it in app folder
compilerstr = Gencode.compilecode(pde);

# run executable file to compute solution and store it in dataout folder
runstr = Gencode.runcode(pde);

# get solution from output files in dataout folder
sol = Postprocessing.fetchsolution(pde,master,dmd);

# visualize the numerical solution of the PDE model using Paraview
pde.visscalars = ["density", 1];  # list of scalar fields for visualization
pde.visvectors = ["velocity", [2, 3]]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");
