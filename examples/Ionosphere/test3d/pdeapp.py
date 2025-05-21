# import external modules
import numpy, os, netCDF4

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
pde['porder'] = 2;         # polynomial degree
pde['torder'] = 2
pde['nstage'] = 2
pde['dt'] = 0.05*numpy.ones(shape=200)
pde['soltime'] = numpy.arange(start=0,stop=200,step=1)
pde['visdt'] = 0.05
pde['pgauss'] = 2*pde['porder'] # number of Gauss points for quadrature
pde['physicsparam'] = numpy.array([1.0, 0.0]);   # unit thermal conductivity
pde['tau'] = numpy.array([1.0]);            # DG stabilization parameter
pde['linearsolvertol'] = 1e-4; # linear solver tolerance
pde['ppdegree'] = 6; # degree of polynomial precontiditoning
pde['RBdim'] = 0; # reduced basis dimension for preconditioner
pde['GMRESrestart'] = 50; # GMRES restart parameter
pde['linearsolveriter'] = 51
pde['precMatrixType'] = 2
pde['NLiter'] = 2

m = 16*1.67e-27
kB = 1.38e-23
pde['mu'] = [m, kB, 1.6e-19, 8.7]

# Choose computing platform and set number of processors
#pde['platform'] = "gpu";   # choose this option if NVIDIA GPUs are available
pde['mpiprocs'] = 1;        # number of MPI processors
pde['hybrid'] = 1;          # 0 -> LDG, 1 -> HDG

Re = 6378
R0 = (Re + 100)/Re
R1 = (Re + 500)/Re
mesh['p'], mesh['t'], mesh['dgnodes'] = Mesh.cubesphere(pde['porder'],R0,R1,15,10)
mesh['boundaryexpr'] = [
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R0**2) < 1e-6,
    lambda p: abs(p[0,:]**2 + p[1,:]**2 + p[2,:]**2 - R1**2) < 1e-6]
mesh['boundarycondition'] = numpy.array([1, 2]); # Set boundary condition for each boundary

n_points_per_element, n_dimensions, n_elements = numpy.shape(mesh["dgnodes"])
vdg = numpy.ndarray(shape=(n_points_per_element, 14, n_elements))
udg = numpy.ndarray(shape=(n_points_per_element, n_dimensions, n_elements))
data = netCDF4.Dataset(filename='input.nc')
vdg[:, 0, :] = data['EX'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 1, :] = data['EY'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 2, :] = data['EZ'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 3, :] = data['BX'][:].reshape((n_points_per_element, n_elements))
vdg[:, 4, :] = data['BY'][:].reshape((n_points_per_element, n_elements))
vdg[:, 5, :] = data['BZ'][:].reshape((n_points_per_element, n_elements))
vdg[:, 6, :] = data['U'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 7, :] = data['V'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 8, :] = data['W'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 9, :] = data['NU'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 10, :] = data['P'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 11, :] = data['L'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 12, :] = data['H'][0, :].reshape((n_points_per_element, n_elements))
vdg[:, 13, :] = data['C'][0, :].reshape((n_points_per_element, n_elements))
udg[:, 0, :] = data['N'][0, :].reshape((n_points_per_element, n_elements))
udg[:, 1, :] = m * (data['N'][0, :] * data['VX'][0, :]).reshape((n_points_per_element, n_elements))
udg[:, 2, :] = m * (data['N'][0, :] * data['VY'][0, :]).reshape((n_points_per_element, n_elements))
udg[:, 3, :] = m * (data['N'][0, :] * data['VZ'][0, :]).reshape((n_points_per_element, n_elements))
udg[:, 4, :] = kB * (data['N'][0, :] * data['T'][0, :]).reshape((n_points_per_element, n_elements))
data.close()

mesh['vdg'] = vdg
mesh['udg'] = udg

# call exasim to generate and run C++ code to solve the PDE model
sol, pde, mesh  = Postprocessing.exasim(pde,mesh)[0:3]

# visualize the numerical solution of the PDE model using Paraview
#pde['paraview'] = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview"; # Set the path to Paraview executable
# pde['visscalars'] = ["temperature", 0]; # list of scalar fields for visualization
# pde['visvectors'] = ["temperature gradient", numpy.array([1, 2, 3]).astype(int)]; # list of vector fields for visualization
# dgnodes = Postprocessing.vis(sol,pde,mesh); # visualize the numerical solution
#x = mesh['dgnodes'][:,0,:]; y = mesh['dgnodes'][:,1,:]; z = mesh['dgnodes'][:,2,:]
#uexact = numpy.sin(numpy.pi*x)*numpy.sin(numpy.pi*y)*numpy.sin(numpy.pi*z); # exact solution
#uh = sol[:,0,:];  # numerical solution
#print("Maximum absolute error: %g\n" % max(abs(uh.flatten()-uexact.flatten())))
print("Done!")
