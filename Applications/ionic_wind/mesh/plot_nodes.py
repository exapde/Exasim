import numpy as np
import matplotlib.pyplot as plt

from read_gmsh import read_gmsh

p,t,pg_nodes,pg_names = read_gmsh('test.msh', 2, True)
print(p.shape)
print(t.shape)
print(pg_nodes.keys())
print(pg_names)

y_coord_nodes_on_symm = np.abs(p[p[:,0]<1e-8][:,1])
y_coord_nodes_on_symm = np.sort(y_coord_nodes_on_symm)

meshsize = np.diff(y_coord_nodes_on_symm)

alpha = 0.5*(np.tanh(100000*(y_coord_nodes_on_symm-0.000391876))+1)
y_predicted = (1-alpha)*(0.047619047619*y_coord_nodes_on_symm+9.52380952e-7) + alpha*(0.074074074074*y_coord_nodes_on_symm- 9.41471029e-6)
y_predicted = np.minimum(y_predicted, np.ones_like(y_predicted)*1e-3)

plt.scatter(y_coord_nodes_on_symm[:-1]*1e3, meshsize*1e3)
plt.plot(y_coord_nodes_on_symm*1e3, y_predicted*1e3, c='red')
plt.axhline(y=0, c='black')
plt.axvline(x=0, c='black')
plt.title('Mesh size progression from needle tip, Chen 2017')
plt.xlabel('Distance from wall [mm]')
plt.ylabel('Element size h [mm]')
plt.show()