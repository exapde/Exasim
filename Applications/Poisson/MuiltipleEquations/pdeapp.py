# Specify an Exasim version to run
version = "Version0.3";

# import external modules
import numpy, os

# Add Exasim to Python search path
cdir = os.getcwd(); ii = cdir.find("Exasim");
exec(open(cdir[0:(ii+6)] + "/Installation/setpath.py").read());

# import internal modules
import Preprocessing, Postprocessing, Gencode, Mesh

# create pde and mesh for each PDE model
pde = [None] * 2
mesh = [None] * 2
exec(open("pdeapp1.py").read());
exec(open("pdeapp2.py").read());

# call exasim to generate and run C++ code to solve the PDE models
sol,pde,mesh,master,dmd,compilerstr,runstr = Postprocessing.exasim(pde,mesh)[0:7];

# visualize the numerical solution of the PDE model using Paraview
for m in range(0,len(pde)):
    pde[m]['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
    pde[m]['visvectors'] = ["temperature gradient", numpy.array([1, 2]).astype(int)]; # list of vector fields for visualization
    pde[m]['visfilename'] = "dataout" + str(m+1) + "/output"; 
    Postprocessing.vis(sol[m],pde[m],mesh[m]); # visualize the numerical solution

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
