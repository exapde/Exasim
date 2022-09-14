import numpy as np
import os


def pdeparams(pde,mesh):

    ##User defined parameters
    #Planet
    planet = "Earth"
    #Set species (or "air" for mixture)
    species = "O"

    #Model (0:Cartesian, 1:cylindrical, 2:spherical)
    model = 2

    #time step (s)
    tstep = 5
    #Length of simulation (days)
    tSimulation = 5
    #frequency of data (minutes)
    freq = 30
    #restart at given time step
    tRestart = 0

    #Coordinates
    longitude = 0
    latitude = 0

    #polynomial order
    porder = 2

    #EUV efficiency
    EUVeff = 1.2

    #Domain of interest
    L = 500e3
    hbot = 100e3
    htop = hbot+L
    lambda0 = 1e-9

    #Initial Conditions
    Tbot = 200
    Ttop = 1000

    ##Read input csv files
    EUV = np.loadtxt('inputs/euv.csv', delimiter=',')
    orbits = np.loadtxt('inputs/neutrals.csv', delimiter=',')
    neutrals = np.loadtxt('inputs/orbits.csv', delimiter=',')



    return pde,mesh

