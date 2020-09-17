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

# Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 2;         # polynomial degree
pde['physicsparam'] = numpy.array([1.0, 0.0, 1.0, 0.0, 0.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([2.0]);            # DG stabilization parameter
pde['linearsolveriter'] = 400;
pde['GMRESrestart'] = 200;
pde['NLtol'] = 1e-4;
pde['GMRESortho'] = 1;
pde['RBdim'] = 1;
pde['NLiter'] = 8;

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 4;        # number of MPI processors

# create a mesh of 8 by 8 by 8 hexes for a unit cube
nd = 3; elemtype = 0;
mesh['p'], mesh['t'] = Mesh.gmshcall(pde, "coneincube", nd, elemtype)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < -10.0+1e-3), lambda p: (p[0,:] > 10.0-1e-3), lambda p: (p[1,:] > 10.0-1e-3), lambda p: (p[0,:] < -10.0+1e-3), lambda p: (p[2,:] < -10.0+1e-3), lambda p: (p[2,:] > 10.0-1e-3), lambda p: (p[2,:] < 1e+3)];
mesh['boundarycondition'] = numpy.array([2, 2, 2, 2, 2, 2, 1]); # Set boundary condition for each boundary
#
# # call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];
#Postprocessing.producecode(pde,mesh);

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["temperature gradient", numpy.array([1, 2, 3]).astype(int)]; # list of vector fields for visualization
dgnodes = Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
# x = dgnodes[:,0,:]; y = dgnodes[:,1,:]; z = dgnodes[:,2,:];
# uexact = numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)*numpy.sin(numpy.pi*z); # exact solution
# uh = sol[:,0,:];  # numerical solution
# print("Maximum absolute error: %g\n" % max(abs(uh.flatten()-uexact.flatten())));
# print("Done!");

# for i in range(0,4):
#     filemesh1 = cdir + "/datain/mesh" + str(i+1) + ".bin";
#     mesh1 = Preprocessing.readmeshstruct(filemesh1);
#     filemesh2 = cdir + "/Applications/Poisson/Cone/datain/mesh" + str(i+1) + ".bin";
#     mesh2 = Preprocessing.readmeshstruct(filemesh2);
#     diff = Preprocessing.checkmesh(mesh1,mesh2);
#     print(diff)
#
#     filemaster1 = cdir + "/facecon" + str(i+1) + ".bin";
#     tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
#     filemaster2 = cdir + "/Applications/Poisson/Cone/facecon" + str(i+1) + ".bin";
#     tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
#     print(max(abs(tm1.flatten('F')+1.0-tm2.flatten('F'))))
#

# filemaster1 = cdir + "/elemrecv1.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/recvelem1.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# n = int(len(tm1)/3.0);
# tm1 = numpy.reshape(tm1,(n,3),'F');
# tm2 = numpy.reshape(tm2,(n,3),'F');
# print(max(abs(tm1.flatten('F')+1.0-tm2.flatten('F'))))
# print(tm1[0:10,:])
# print(tm2[0:10,:])
#
# filemaster1 = cdir + "/elemsend1.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/elemsend1.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')+1.0-tm2.flatten('F'))))
# print(tm1[0:10])
# print(tm2[0:10])

# filemaster1 = cdir + "/t.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/t.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')+1.0-tm2.flatten('F'))))
#
# filemaster1 = cdir + "/t2t.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/t2t.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')+1.0-tm2.flatten('F'))))
# print(tm1.shape)
# print(tm2.shape)
# print(tm1[0:10])
# print(tm2[0:10])

# fileapp1 = cdir + "/datain/app.bin";
# app1 = Preprocessing.readapp(fileapp1);
# fileapp2 = cdir + "/Applications/Poisson/Cone/datain/app.bin";
# app2 = Preprocessing.readapp(fileapp2);
# diff = Preprocessing.checkapp(app1,app2);
# print(app1['problem'])
# print(app2['problem'])
# print(diff)
#
# filemaster1 = cdir + "/datain/master.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/datain/master.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
#
# filemesh1 = cdir + "/datain/mesh1.bin";
# mesh1 = Preprocessing.readmeshstruct(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Cone/datain/mesh1.bin";
# mesh2 = Preprocessing.readmeshstruct(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)
#
# filemesh1 = cdir + "/datain/mesh2.bin";
# mesh1 = Preprocessing.readmeshstruct(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Cone/datain/mesh2.bin";
# mesh2 = Preprocessing.readmeshstruct(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)
#
# filemesh1 = cdir + "/datain/mesh3.bin";
# mesh1 = Preprocessing.readmeshstruct(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Cone/datain/mesh3.bin";
# mesh2 = Preprocessing.readmeshstruct(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)
#
# filemesh1 = cdir + "/datain/mesh4.bin";
# mesh1 = Preprocessing.readmeshstruct(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Cone/datain/mesh4.bin";
# mesh2 = Preprocessing.readmeshstruct(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)


# filemaster1 = cdir + "/datain/sol1.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/datain/sol1.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
#
# filemaster1 = cdir + "/datain/sol2.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/datain/sol2.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
#
# filemaster1 = cdir + "/datain/sol3.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/datain/sol3.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
#
# filemaster1 = cdir + "/datain/sol4.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Cone/datain/sol4.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))
