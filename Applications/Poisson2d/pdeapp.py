# Specify an Exasim version to run
version = "Version0.1";

# import external modules
import numpy, os, sys

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
versiondir = cdir[0:(ii+6)] + "/"  + version + "/Python";
exec(open(versiondir + "/setpath.py").read());
# n = len(sys.path);
# for i in range(0,n):
#     print(sys.path[i])

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# Create app struct
app = Preprocessing.initializeapp(version);

# Define a PDE model: governing equations and boundary conditions
app['appname'] = "poisson";       # application name
app['pdemodel'] = "ModelD";       # ModelC, ModelD, ModelW
app['pdemodelfile'] = "pdemodel"; # name of a file defining the PDE model

# Set compilers and compiling options
app['platform'] = "cpu";    # choose this option if NVIDIA GPUs are not available
#app['platform = "gpu";   # choose this option if NVIDIA GPUs are available
app['mpiprocs'] = 2;        # number of MPI processors
app['cpucompiler'] = "g++"; # Clang/GNU C++ compiler
if app['mpiprocs']>1:        # set MPI compiler if app['mpiprocs > 1
    app['mpicompiler'] = "/opt/local/bin/mpicxx-openmpi-mp"; # path to a MPI compiler
    app['mpirun'] = "/opt/local/bin/mpirun-openmpi-mp";
if app['platform'] == "gpu":
    app['gpucompiler'] = "nvcc"; # Nvidia CUDA compiler (AMD GPU will be supported in the next version)
Gencode.checkcompilers(app);     # check if the compilers are available

# Set discretization parameters, physical parameters, and solver parameters
app['porder'] = 3;         # polynomial degree
app['dt'] = numpy.array([0.0]);             # steady-state problem
app['boundaryconditions'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary
app['physicsparam'] = numpy.array([1.0]);   # unit thermal conductivity
app['tau'] = numpy.array([1.0]);            # DG stabilization parameter
app['ncu'] = 1;

# create a linear mesh
m = 8; n = 8; # m by n grid of the square
elemtype = 1; # quad elements
mesh = {'p' : [], 't' : []};
mesh['p'], mesh['t'] = Mesh.SquareMesh(m,n,elemtype)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3)];
# expressions for curved boundaries
mesh['curvedboundary'] = [0, 0, 0, 0];
mesh['curvedboundaryexpr'] = [];
# experssions for periodic boundaries
mesh['periodicexpr'] = [];

# generate input files and store them in datain folder
app, mesh, master, dmd = Preprocessing.preprocessing(app,mesh);

# generate source codes and store them in app folder
Gencode.gencode(app);

# compile source codes to build an executable file and store it in app folder
compilerstr = Gencode.compilecode(app);

# run executable file to compute solution and store it in dataout folder
runstr = Gencode.runcode(app);

# get solution from output files in dataout folder
app['vistime'] = [];
sol = Postprocessing.fetchsolution(app,master,dmd);

# perform visualization
app['paraview'] = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview";
app['visfilename'] = "dataout/output";
app['visscalars'] = ["temperature", 0];
app['visvectors'] = ["temperature gradient", numpy.array([1, 2]).astype(int)];
app['viselem'] = numpy.arange(0,mesh['t'].shape[1]);
Postprocessing.vis(sol,app,mesh,master); # visualize the numerical solution
print("Done!");





















# fileapp1 = cdir + "/datain/app['bin";
# app1 = Preprocessing.readapp(fileapp1);
# fileapp2 = cdir + "/Applications/Poisson2d/datain/app['bin";
# app2 = Preprocessing.readapp(fileapp2);
# diff = Preprocessing.checkapp(app1,app2);
# print(diff)
#
# filemaster1 = cdir + "/datain/master.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson2d/datain/master.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
#
# filemesh1 = cdir + "/datain/mesh.bin";
# mesh1 = Preprocessing.readmesh(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson2d/datain/mesh.bin";
# mesh2 = Preprocessing.readmesh(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(diff)

# tm1 = numpy.fromfile(open(filemesh1, "r"), dtype=numpy.float64);
# tm2 = numpy.fromfile(open(filemesh2, "r"), dtype=numpy.float64);
# print(mesh1['nsize'])
# print(mesh2['nsize'])
# k1 = 0; k2 = 20;
# print(max(abs(tm1[k1:k2].flatten('F')-tm2[k1:k2].flatten('F'))))
# k1 = 20; k2 = 1152+20;
# print(max(abs(tm1[k1:k2].flatten('F')-tm2[k1:k2].flatten('F'))))
# print(tm1[k1:k2])
# print(tm2[k1:k2])
# print(mesh1['facecon'].flatten('F'))
# print(mesh2['facecon'].flatten('F'))

# print(tm1.shape)
# print(tm2.shape)

# print(mesh1['colent2elem'].T)
# print(mesh2['colent2elem'].T)

# print(mesh['f'].T)
# print(mesh['dgnodes'][:,:,0])
# print(mesh['dgnodes'][:,:,-1])
