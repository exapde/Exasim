function [pde, mesh] = pdeparams(pde,mesh)

%% User defined parameters
%Planet
planet = "Earth";
%Set species (or "air", for mixture)
species = "O";

%Model (0:Cartesian, 1:cylindrical, 2:spherical)
model = 2;

%timeStep (s)
tstep = 5;
%Length of simulation (days)
tSimulation = 5;
%frequency of data (minutes)
freq = 30;
%Restart at given time step
tRestart = 0;

%Coordinates
longitude = 0;
latitude = 0;

%Polynomial order
porder = 2;

%EUV efficiency
EUVeff = 1.2;

%Domain of interest
L = 500e3;
hbot = 100e3;     %height inner surface (km)
htop = hbot+L;     %height outer surface (km)
lambda0 = 1e-9;   %reference wavelength (1nm=1e-9m)

%Initial condition
Tbot = 200;     %(in K)
Ttop = 1000;    %(in K)


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
speciesEUV = string(table2array(EUV(5:end,2)));
iSpecies = find(strcmp(string(table2array(neutrals(:,1))),species));
iSpeciesEUV = find(strcmp(speciesEUV,species));
neutralSpecies = string(table2array(neutrals(:,1)));
neutrals = table2array(neutrals(:,2:end));
gam = 5/3;
amu = 1.66e-27;
if strcmp(species,"air")
    rho = (neutrals(:,1)'*neutrals(:,end))*amu;
    m = rho/sum(neutrals(:,end));
    expKappa = (neutrals(:,4)'*neutrals(:,end))/sum(neutrals(:,end));
    kappa0 = (neutrals(:,3)'*neutrals(:,end))/sum(neutrals(:,end))*(Tbot^expKappa);
    
    iEUV = zeros(size(speciesEUV));
    for iSpecies=1:length(speciesEUV)
        index = find(strcmp(speciesEUV(iSpecies),neutralSpecies));
        if isempty(index)
            iEUV(iSpecies) = 0;
        else
            iEUV(iSpecies) = index;
        end
    end
    localIndicesEUV = find(iEUV);
    crossSections_d = neutrals(iEUV(localIndicesEUV),end)'*(table2array(EUV(iEUV(localIndicesEUV)+4,6:42)).*table2array(EUV(iEUV(localIndicesEUV)+4,4)))/sum(neutrals(:,end));     %m2
else
    m = neutrals(iSpecies,1)*amu;
    rho = neutrals(iSpecies,end)*m;
    expKappa = neutrals(iSpecies,4);
    kappa0 = neutrals(iSpecies,3)*(Tbot^expKappa);
    crossSections_d = table2array(EUV(iSpeciesEUV+4,6:42))*table2array(EUV(iSpeciesEUV+4,4));     %m2
end
lambda_d = (0.5*(table2array(EUV(1,6:42))+table2array(EUV(2,6:42))))*1e-10;    % initially in Armstrongs
AFAC = table2array(EUV(4,6:42));
F74113_d = table2array(EUV(3,6:42))*table2array(EUV(3,4))*1e4;    %1/(m2*s)
clear neutrals
clear EUV

%% Physical quantities
kBoltzmann = 1.38e-23;
hPlanck = 6.0626e-34;
c = 3e8;
gravitationalConstant = 6.674e-11;
R = kBoltzmann/m;
g = gravitationalConstant*planetMass/radiusIn^2;
omega = 2*pi/periodDay;
cp = gam*R/(gam-1);
H = R*Tbot/g;

%Reference quantities
T0 = 1;
T1 = Ttop/Tbot;
epsilon = abs(T1-T0);
v0 = sqrt(gam*R*Tbot);
t0 = H/v0;
R0 = radiusIn/H;
R1 = radiusOut/H;

expMu = 0.5;
expKappa = 0.75;

alpha0 = kappa0/(rho*cp);
mu0 = 1.3e-4*(Tbot/R)^expMu;
nu0 = mu0/rho;

lambda = lambda_d/lambda0;
crossSections = crossSections_d/H^2;
F74113 = F74113_d*(H^2*t0);

tauA = 5;

%Nondimensional quantities
Gr = g*H^3/nu0^2;
Pr = nu0/alpha0;
Fr = sqrt(omega^2*H/g);

Keuv = (gam*kBoltzmann*Tbot)/(hPlanck*c/lambda0);
M = rho*H^3/m;


%% Time parameters;
tstepStar = tstep/t0;
nTimeSteps = ceil(tSimulation*periodDay/tstep);
freqTimeSteps = ceil(freq*60/tstep);

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = porder;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages
pde.dt = tstepStar*ones(360,1);   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = freqTimeSteps;          % solution is saved every 100 time steps
pde.soltime = freqTimeSteps:freqTimeSteps:length(pde.dt); % steps at which solution are collected
pde.timestepOffset = tRestart;


%% Vectors of physical and external parameters
pde.physicsparam = [gam Gr Pr Fr Keuv  M  rho T0 T1 R0 R1 H EUVeff model longitude latitude declinationSun tauA t0];
                   % 1  2  3  4   5    6   7  8  9  10 11 12 13     14      15      16           17         18  19

%External params (EUV)
pde.externalparam = [lambda,crossSections,AFAC,F74113];


%% Solver parameters
pde.extStab = 1;
pde.tau = 0.0;                  % DG stabilization parameter
pde.GMRESrestart=29;            % number of GMRES restarts
pde.linearsolvertol=1e-16;     % GMRES tolerance
pde.linearsolveriter=30;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-10;               % Newton toleranccd dataoue
pde.NLiter = 2;                 % Newton iterations
pde.matvectol = 1e-7;
pde.RBdim = 8;


%% Grid
[mesh.p,mesh.t] = mesh1D_adapted(R0,R1);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:)-R0)<1e-8, @(p) abs(p(1,:)-R1)<1e-8};
mesh.boundarycondition = [1;2];  % Left, Right

