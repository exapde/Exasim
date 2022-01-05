from initializepde import initializepde
from initializemesh import initializemesh

def initializeexasim():

    version = "src";

    pde = initializepde(version);

    mesh = initializemesh(version);

    return pde, mesh
