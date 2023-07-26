function [pde, mesh] = pdeparams(pde,mesh)
addpath('utilitiesMSIS/');

%Example and case
example = 2;
section = 1;

%Mesh refinement
iH = 3;
%Polynomial order
porder = 2;

switch porder
    case 1
        tstep = 5;
    case 2
        tstep = [2.5,5];
    otherwise
        tstep = [2,4];    
end

%Problem setup
[date,tSimulation,F10p7,F10p7a,Ap,L] = setupCases(example,section);

%% User defined parameters
%Planet
planet = "Earth";
%Set species 
nspecies = 4;

%start of simulation
% date = [2002, 3, 18];       %% year month day
hour = [00, 00, 00];        %% HH MM SS (24h format)

%Length of simulation (days)
% tSimulation = 5;
%timeStep (s)
% tstep = [5,10];
tpoints = [0];    %in hours
%frequency of data (minutes)
freq = tstep/60;
%Restart at given time step
% tRestart = 31680;
tRestart = 0;

%Stabilization parameter
tauA = 10;

%EUV parameters
alpha0 = 1;
EUVeff = 0.25;
% F10p7 = 150;
% F10p7a = 150;

%Domain of interest
% L = 500e3;
hbot = 100e3;     %height inner surface (km)
htop = hbot+L;     %height outer surface (km)
lambda0 = 1e-9;   %reference wavelength (1nm=1e-9m)

%Calculations
year = date(1);
doy = day(datetime(date),'dayofyear');
sec = hour(3) + 60*(hour(2) + 60*hour(1));
Nday = doy + sec/86400;

%Define species
species = ["O"; "N2"; "O2"; "He"];
indicesMSIS = [2;3;4;1];
species = species(1:nspecies);
indicesMSIS = indicesMSIS(1:nspecies);

%% Read input csv files
addpath('inputs');
load('inputCSV.mat','orbits','neutrals','EUV');

%Planet information
iPlanet = strcmp(string(table2array(orbits(:,1))),planet);
periodDay = table2array(orbits(iPlanet,14))*3600;
radius = table2array(orbits(iPlanet,18))*1e3;
radiusIn = radius+hbot;
radiusOut= radius+htop;
planetMass = table2array(orbits(iPlanet,17));
declinationSun = table2array(orbits(iPlanet,20));
clear orbits

%Species information
iSpecies = zeros(nspecies,1);
iSpeciesEUV = zeros(nspecies,1);
for isp = 1:nspecies
    iSpecies(isp) = find(strcmp(string(table2array(neutrals(:,1))),species(isp)));
    iSpeciesEUV(isp) = find(strcmp(table2array(EUV(:,2)),species(isp)));
end

amu = 1.66e-27;
mass = table2array(neutrals(iSpecies,2))*amu;
ckappa0 = table2array(neutrals(iSpecies,4));% + 2e-4;
% expKappa = table2array(neutrals(iSpecies(1),5));
expKappa = 0.69;

lambda_d = (0.5*(table2array(EUV(1,6:42))+table2array(EUV(2,6:42))))*1e-10;    % initially in Armstrongs
AFAC = table2array(EUV(4,6:42));
F74113_d = table2array(EUV(3,6:42))*table2array(EUV(3,4))*1e4;    %1/(m2*s)
crossSections_d = table2array(EUV(iSpeciesEUV,6:42)).*table2array(EUV(iSpeciesEUV,4));     %m2

clear neutrals
clear EUV

%% Reference values
m = mass(1);
gam = 5/3;
%MSIS reference values
[rho0,T0,chi,cchi] = MSIS_referenceValues(indicesMSIS,mass,hbot,htop,year,doy,sec,F10p7,F10p7a,Ap);

%% Physical quantities
kBoltzmann = 1.38e-23;
hPlanck = 6.0626e-34;
c = 3e8;
gravitationalConstant = 6.674e-11;
R = kBoltzmann/m;
g = gravitationalConstant*planetMass/radiusIn^2;
omega = 2*pi/periodDay;
cp = gam*R/(gam-1);
H = R*T0/g;

%Reference quantities
v0 = sqrt(gam*R*T0);
t0 = H/v0;
R0 = radiusIn/H;
R1 = radiusOut/H;

%Viscosity and conductivity
cmu0 = 5*1.3e-4;
expMu = 0.5;
ckappai = ckappa0/(chi(1,:)*ckappa0);
kappa0 = chi(1,:)*ckappa0*T0^expKappa;
alpha0 = alpha0*kappa0/(rho0*cp);
mu0 = cmu0*(T0/R)^expMu;
nu0 = mu0/rho0;

nuEddy = 200;
alphaEddy = 10;

%Parameters for EUV
lambda = lambda_d/lambda0;
crossSections = crossSections_d/H^2;
F74113 = F74113_d*(H^2*t0);
mass = mass/m;

%Nondimensional quantities
Gr = g*H^3/nu0^2;
Pr = nu0/alpha0;
Fr = sqrt(omega^2*H/g);

Keuv = (gam*kBoltzmann*T0)/(hPlanck*c/lambda0);
M = rho0*H^3/m;

%% Time parameters;
if tRestart>0
    t0real = (tRestart-3600*tpoints(1)/tstep(1))*tstep(2) + 3600*tpoints(1);
    tstep = tstep(2);
    tpoints = [];
else
    t0real = 0;
end
tSimulation = tSimulation - t0real/86400;

tend = tSimulation*periodDay;
tpoints = unique([0,tpoints*3600,tend]);
tintervals = tpoints(2:end)-tpoints(1:end-1);
nTimeStepsI = ceil(tintervals./tstep);
dt = [];
soltime = [];
soltime0 = 0;
freqTimeSteps = ceil(freq*60./tstep);
for it = 1:length(tintervals)
    dt = [dt; tstep(it)/t0*ones(nTimeStepsI(it),1)];
    soltime = [soltime, soltime0 + freqTimeSteps(it):freqTimeSteps(it):nTimeStepsI(it)];
    soltime0 = soltime(end);
end

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = porder;          % polynomial degree
pde.torder = max(porder,2);          % time-stepping order of accuracy
pde.nstage = max(porder,2);          % time-stepping number of stages
pde.dt = dt;   % time step sizes
pde.visdt = freq*60;         % visualization timestep size
pde.saveSolFreq = min(freqTimeSteps);          % solution is saved every 100 time steps
pde.soltime = soltime; % steps at which solution are collected
pde.time = (t0real - tRestart*tstep(1))/t0;
pde.timestepOffset = tRestart;


%% Vectors of physical and external parameters
pde.physicsparam = [gam Gr Pr Fr Keuv  M EUVeff declinationSun F10p7 F10p7a Nday expMu expKappa nuEddy alphaEddy R0 R1 H T0 rho0 t0, tauA, year];
                   
% External params (EUV)
pde.externalparam = [lambda,AFAC,F74113,reshape(crossSections',[37*nspecies,1])',reshape(cchi',[4*(nspecies-1),1])',mass',ckappai'];        

% Parameters for MSIS calls
paramsMSIS = [R0,R1,year,doy,sec,F10p7,F10p7a,hbot,H,T0,rho0,Fr,m,Ap];


%% Solver parameters
pde.extStab = 1;
pde.tau = 0.0;                  % DG stabilization parameter
pde.GMRESrestart=19;            % number of GMRES restarts
pde.linearsolvertol=1e-16;      % GMRES tolerance
pde.linearsolveriter=20;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-10;              % Newton tolerance
pde.NLiter = 4;                 % Newton iterations
pde.matvectol = 1e-4;
pde.RBdim = 10;
pde.matvecorder = 2;

%% Grid
[mesh.p,mesh.t,mesh.dgnodes] = sphereCubeRefinement(pde.porder, R0, R1, iH);
% [mesh.p,mesh.t, mesh.dgnodes] = sphereCubeAdaptedScaleHeight(pde.porder, R0, R1,paramsMSIS,indicesMSIS,mass);
% [mesh.p,mesh.t, mesh.dgnodes] = sphereCubeAdapted(pde.porder, R0, R1);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:).^2+p(2,:).^2+p(3,:).^2-R0^2)<1e-6, @(p) abs(p(1,:).^2+p(2,:).^2+p(3,:).^2-R1^2)<1e-6};
mesh.boundarycondition = [1 2];  % Inner, Outer

if tRestart<1
    %% Initial condition
    [npoint,ndim,nelem] = size(mesh.dgnodes);
    ndg = npoint*nelem;
    ncu = 5;
    nc = ncu*(ndim+1);
    xdg = reshape(pagetranspose(mesh.dgnodes),[ndim,ndg])';

%     u0 = MSIS_initialCondition3D(xdg,paramsMSIS,indicesMSIS,mass);
    u0 = MSIS_hydrostaticInitialCondition3D(xdg,paramsMSIS,indicesMSIS,mass);
    u0 = pagetranspose(reshape(u0',[nc,npoint,nelem]));
    mesh.udg = u0;
end

