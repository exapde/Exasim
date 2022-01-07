function initializeexasim(version)

# create pde model
pde = Preprocessing.initializepde(version);

# create a mesh structure
mesh = Preprocessing.initializemesh(version);

return pde, mesh

end
