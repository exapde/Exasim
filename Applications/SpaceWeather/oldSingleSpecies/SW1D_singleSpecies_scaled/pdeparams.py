from numpy import *
from pandas import *
from math import *
import os
from mesh1D_adapted import mesh1D_adapted 

def pdeparams(pde,mesh):

    ##User defined parameters
    #Planet
    planet = 'Earth'
    #Set species (or "air" for mixture)
    species = 'O'

    #Model (0:Cartesian, 1:cylindrical, 2:spherical)
    model = 0

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
    EUV = read_csv('inputs/euv.csv',header=None)
    orbits = read_csv('inputs/orbits.csv')
    neutrals = read_csv('inputs/neutrals.csv',delimiter=";")

    #Planet information
    iPlanet = orbits.values[:,0]==planet
    periodDay = float(orbits.values[iPlanet,13])*3600
    radius = float(orbits.values[iPlanet,17])*1e3
    radiusIn = radius+hbot
    radiusOut = radius+htop
    planetMass = float(orbits.values[iPlanet,16])
    declinationSun = float(orbits.values[iPlanet,19])

    #Species information
    speciesEUV = EUV.values[4:,1]
    iSpecies = where(neutrals.values[:,0]==species)[0]
    iSpeciesEUV = where(speciesEUV==species)[0]
    neutralSpecies = neutrals.values[:,0]
    neutrals = neutrals.values[:,1:]
    gam = 5/3
    amu = 1.66e-27
    if species=='air':
        m = 1
        rho = 2
    else:
        m = neutrals[iSpecies,0][0]*amu
        rho = neutrals[iSpecies,-1][0]*m
        expKappa = neutrals[iSpecies,3][0]
        kappa0 = neutrals[iSpecies,2][0]*(Tbot**expKappa)
        crossSections_d = EUV.values[iSpeciesEUV+4,5:42]*float(EUV.values[iSpeciesEUV+4,3])

    lambda_d = 0.5*(EUV.values[0,5:42]+EUV.values[1,5:42])*1e-10
    AFAC = EUV.values[3,5:42]
    F74113_d = EUV.values[2,5:42]*float(EUV.values[2,3])*1e4

    ##Physical quantities
    kBoltzmann = 1.38e-23
    hPlanck = 6.0626e-34
    c = 3e8
    gravitationalConstant = 6.674e-11
    R = kBoltzmann/m
    g = gravitationalConstant*planetMass/radiusIn**2
    omega = 2*pi/periodDay
    cp = gam*R/(gam-1)
    H = R*Tbot/g

    ##Reference quantities
    T0 = 1
    T1 = Ttop/Tbot
    epsilon = abs(T1-T0)
    v0 = sqrt(gam*R*Tbot)
    t0 = H/v0
    R0 = radiusIn/H
    R1 = radiusOut/H

    expMu = 0.5
    expKappa = 0.75

    alpha0 = kappa0/(rho*cp)
    mu0 = 1.3e-4*(Tbot/R)**expMu
    nu0 = mu0/rho

    lambda_EUV = lambda_d/lambda0
    crossSections = crossSections_d/H**2
    F74113 = F74113_d*(H**2*t0)

    tauA = 10.0

    ## Nondimensional quantities
    Gr = g*H**3/nu0**2
    Pr = nu0/alpha0
    Fr = sqrt(omega**2*H/g)

    Keuv = (gam*kBoltzmann*Tbot)/(hPlanck*c/lambda0)
    M = rho*H**3/m

    ## Time parameters
    tstepStar = tstep/t0
    nTimeSteps = ceil(tSimulation*periodDay/tstep)
    freqTimeSteps = ceil(freq*60/tstep)

    # Set discretization parameters, physical parameters, and solver parameters
    pde['porder'] = porder                     #polynomial degree
    pde['torder'] = 2                           #time-stepping order of accuracy
    pde['nstage'] = 2                          #time-stepping number of stages
    pde['dt'] = tstepStar*ones([nTimeSteps,1])   #time step sizes
    pde['visdt'] = pde['dt'][0]                   #visualization timestep size
    pde['saveSolFreq'] = freqTimeSteps         #solution is saved every 100 time steps
    pde['soltime'] = arange(freqTimeSteps,pde['dt'].shape[0],freqTimeSteps) #steps at which solution are collected
    pde['timestepOffset'] = tRestart


    pde['physicsparam'] = array([gam,Gr,Pr,Fr,Keuv,M,rho,T0,T1,R0,R1,H,EUVeff,model,longitude,latitude,declinationSun,tauA,t0])
    pde['externalparam'] = hstack([lambda_EUV,crossSections[0,:],AFAC,F74113])

    # Solver parameters
    pde['extStab'] = 1
    pde['tau'] = 0.0  # DG stabilization parameter
    pde['GMRESrestart'] = 29  # number of GMRES restarts
    pde['linearsolvertol'] = 1e-16  # GMRES tolerance
    pde['linearsolveriter'] = 30  # number of GMRES iterations
    pde['precMatrixType'] = 2  # preconditioning type
    pde['NLtol'] = 1e-10  # Newton toleranccd dataoue
    pde['NLiter'] = 2  # Newton iterations
    pde['matvectol'] = 1e-7
    pde['RBdim'] = 8

    # Mesh
    resolution = 16
    mesh['p'], mesh['t'] = mesh1D_adapted(R0,R1,resolution)
    # expressions for domain boundaries
    mesh['boundaryexpr'] = [lambda p: (p[0,:] < R0+1e-3), lambda p: (p[0,:] > R1-1e-3)]
    mesh['boundarycondition'] = array([1, 2])  # Set boundary condition for each boundary

    return pde,mesh

