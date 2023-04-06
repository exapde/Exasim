# Specify an Exasim version to run
version = "Version0.3";

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
pde['mpiprocs'] = 1;        # number of MPI processors

# # Set discretization parameters, physical parameters, and solver parameters
pde['porder'] = 3;         # polynomial degree
pde['physicsparam'] = numpy.array([.01, 0, 1.0]);
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter

pde['torder'] = 2;          # time-stepping order of accuracy
pde['nstage'] = 2;          # time-stepping number of stages
pde['dt'] = 0.02*numpy.ones(100);   # time step sizes
pde['soltime'] = numpy.arange(1,pde['dt'].size+1,1); # steps at which solution are collected
pde['visdt'] = 1.0; # visualization timestep size

pde['GMRESrestart']=15;            # number of GMRES restarts
pde['linearsolvertol']=1e-12;      # GMRES tolerance
pde['linearsolveriter']=16;        # number of GMRES iterations
pde['precMatrixType']=2;           # preconditioning type
pde['NLtol'] = 1e-12;              # Newton tolerance
pde['NLiter']=2;                   # Newton iterations


# # create a mesh of 8 by 8 quads on a square domain
# mesh['p'], mesh['t'] = Mesh.SquareMesh(16,16,1)[0:2];
# # expressions for domain boundaries
# mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3)];
# mesh['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary

# create a mesh of 8 by 8 quads on a square domain
mesh['p'], mesh['t'] = Mesh.SquareMesh(8,8,1)[0:2];
# expressions for domain boundaries
mesh['boundaryexpr'] = [lambda p: (p[1,:] < 1e-3), lambda p: (p[0,:] > 1-1e-3), lambda p: (p[1,:] > 1-1e-3), lambda p: (p[0,:] < 1e-3)];
mesh['boundarycondition'] = numpy.array([1, 1, 1, 1]); # Set boundary condition for each boundary
# mesh['periodicexpr'] = [[2, lambda p: p[1,:], 4, lambda p: p[1,:]], [1, lambda p: p[0,:], 3, lambda p: p[0,:]]];
mesh['periodicexpr'] = [[1, lambda p: p[0,:], 3, lambda p: p[0,:]]];



# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3];

# visualize the numerical solution of the PDE model using Paraview
pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
pde['visvectors'] = ["temperature gradient", numpy.array([1, 2]).astype(int)]; # list of vector fields for visualization
Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
print("Done!");













































# npf = dmd[0]['facecon'].shape[0];
# nf = dmd[0]['facecon'].shape[2];
# print(numpy.reshape(dmd[0]['facecon'][:,0,:],(npf,nf),'F').T)
# print(numpy.reshape(dmd[0]['facecon'][:,1,:],(npf,nf),'F').T)

# fileapp1 = cdir + "/datain/app.bin";
# app1 = Preprocessing.readapp(fileapp1);
# fileapp2 = cdir + "/Applications/Poisson/Poisson2d/datain/app.bin";
# app2 = Preprocessing.readapp(fileapp2);
# diff = Preprocessing.checkapp(app1,app2);
# print(app1['problem'])
# print(app2['problem'])
# print(diff)

# filemaster1 = cdir + "/datain/master.bin";
# tm1 = numpy.fromfile(open(filemaster1, "r"), dtype=numpy.float64);
# filemaster2 = cdir + "/Applications/Poisson/Poisson2d/datain/master.bin";
# tm2 = numpy.fromfile(open(filemaster2, "r"), dtype=numpy.float64);
# print(max(abs(tm1.flatten('F')-tm2.flatten('F'))))

# filemesh1 = cdir + "/datain/mesh1.bin";
# mesh1 = Preprocessing.readmesh(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Poisson2d/datain/mesh1.bin";
# mesh2 = Preprocessing.readmesh(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)
#
# filemesh1 = cdir + "/datain/mesh2.bin";
# mesh1 = Preprocessing.readmesh(filemesh1);
# filemesh2 = cdir + "/Applications/Poisson/Poisson2d/datain/mesh2.bin";
# mesh2 = Preprocessing.readmesh(filemesh2);
# diff = Preprocessing.checkmesh(mesh1,mesh2);
# print(mesh1['nbsd'])
# print(mesh2['nbsd'])
# print(diff)

# print(mesh1['ndims'])
# print(mesh2['ndims'])
# print(mesh1['nsize'])
# print(mesh2['nsize'])
# print(mesh1['facecon'][0:10])
# print(mesh2['facecon'][0:10])

# print(dmd[0]['facecon'][:,:,0])






















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
