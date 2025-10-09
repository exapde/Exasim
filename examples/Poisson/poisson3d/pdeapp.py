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
pde['model'] = "ModelD";       # ModelC, ModelD, ModelW
pde['modelfile'] = "pdemodel"; # name of a file defining the PDE model

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['pgauss'] = 2*pde['porder'] # number of Gauss points for quadrature
pde['physicsparam'] = numpy.array([1.0, 0.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter
pde['linearsolvertol'] = 1e-8; # linear solver tolerance
pde['ppdegree'] = 1; # degree of polynomial precontiditoning
pde['RBdim'] = 0; # reduced basis dimension for preconditioner
pde['GMRESrestart'] = 50; # GMRES restart parameter
pde['linearsolveriter'] = 1000; 
pde['preconditioner'] = 1; 

# Choose computing platform and set number of processors
pde['platform'] = "gpu"
pde['mpiprocs'] = 4
pde['cpucompiler']  = "CC"
pde['mpicompiler'] = "CC"
pde['gpucompiler'] = cdir[0:(ii+6)] + "/kokkos/bin/nvcc_wrapper" # e.g. /path_to_kokkos/bin/nvcc_wrapper

# create a mesh of 8 by 8 by 8 hexes for a unit cube
mesh['p'], mesh['t'] = Mesh.cubemesh(16,16,16,1)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3), lambda p: (p[2,:] < 1e-3), lambda p: (p[2,:] > 1-1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1, 1, 1]); # Set boundary condition for each boundary

# call exasim to generate and run C++ code to solve the PDE model
# sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# search compilers and set options
pde = Gencode.setcompilers(pde)

# generate input files and store them in datain folder
pde, mesh, master, dmd = Preprocessing.preprocessing(pde,mesh)

# generate source codes and store them in app folder
Gencode.gencode(pde)

# compile source codes to build an executable file and store it in build folder
compilerstr = Gencode.cmakecompile(pde)

# Run code
pde['mpirun'] = "mpirun"; # command to run MPI programs
runstr = Gencode.runcode(pde, 1)

pde['vistime'] = [];
sol = Postprocessing.fetchsolution(pde,master,dmd, pde['buildpath'] + "/dataout");

# visualize the numerical solution of the PDE model using Paraview
#pde['paraview'] = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview"; # Set the path to Paraview executable
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["temperature gradient", numpy.array([1, 2, 3]).astype(int)]; # list of vector fields for visualization
dgnodes = Preprocessing.createdgnodes(mesh['p'],mesh['t'],mesh['f'],mesh['curvedboundary'],mesh['curvedboundaryexpr'],pde['porder']);
x = dgnodes[:,0,:]; y = dgnodes[:,1,:]; z = dgnodes[:,2,:];
uexact = numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)*numpy.sin(numpy.pi*z); # exact solution
uh = sol[:,0,:];  # numerical solution
print("Maximum absolute error: %g\n" % max(abs(uh.flatten()-uexact.flatten())));
print("Done!");
