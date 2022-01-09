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

