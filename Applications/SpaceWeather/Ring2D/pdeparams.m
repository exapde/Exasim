function [pde, mesh] = pdeparams(pde,mesh)

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model


% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;          % polynomial degree
pde.torder = 2;          % time-stepping order of accuracy
pde.nstage = 2;          % time-stepping number of stages

dt0 = 5;       %seconds
tend = 1*86400;       %in s
freqT = 5*60;      %in s
pde.dt = dt0*ones(tend/dt0,1);   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = freqT/dt0;          % solution is saved every X time steps
pde.soltime = freqT/dt0:freqT/dt0:length(pde.dt); % steps at which solution are collected
pde.timestepOffset = 0;

R_Earth = 6370;
H0 = 100;
H1 = 600;
R0 = (R_Earth+H0);
R1 = (R_Earth+H1);

gam = 1.4;                      % specific heat ratio
gam1 = gam - 1;   
Minf = 0.05;                     % Mach number (not used actually)

rbot = 660;                    % kg/km3
Tbot = 200;                    % K
Ttop = 1000;                   % K

cp = 1000*1e-6;                 % km2/(s2*K)
pbot = rbot*Tbot*gam1*cp/gam;   % kg/(km*s2)

gravity = 9.5*1e-3;        % km/s^2
omega = 2*pi/86400;        % rad/s

visc0 = 1.3e-4;         % kg/(s*K)
kappa0 = 5.6e-7;        % (*T^0.75) kg*km/(s3*K)
m = 28.9*1.66e-27;      % kg
h = 6.0626e-34*1e-6;    % kg*km2/s   
c = 3e8*1e-3;           % km/s

EUV = readtable('euv.csv');
lambda = (0.5*(table2array(EUV(1,6:42))+table2array(EUV(2,6:42))))*1e-13;                   % km
crossSections = (0.78*table2array(EUV(5,6:42))+0.22*table2array(EUV(8,6:42)))*1e-22*1e-6;   % km2
AFAC = table2array(EUV(4,6:42));                                                            % non-dimensional
F74113 = table2array(EUV(3,6:42))*1e9*1e10;                                                 % 1/(km2*s)
pde.externalparam = [lambda,crossSections,AFAC,F74113];

pde.physicsparam = [gam visc0 kappa0 Minf gravity omega cp rbot pbot Tbot Ttop R0 R1 m h c];
                   % 1    2     3     4      5      6   7   8    9    10   11  12 13 14 15 16 

pde.tau = 10;                  % DG stabilization parameter
pde.GMRESrestart=25;            % number of GMRES restarts
pde.linearsolvertol=1e-3;     % GMRES tolerance
pde.linearsolveriter=100;        % number of GMRES iterations
pde.precMatrixType=2;           % preconditioning type
pde.NLtol = 1e-3;               % Newton tolerance
pde.NLiter = 5;                 % Newton iterations


% Generate grid
rate = 1/3;
AR = 5;
% [mesh.p,mesh.t, mesh.dgnodes] = mkmesh_ring_scaleheight(pde.porder,R0,R1,Tbot,Ttop,gravity,omega,(gam-1)*cp/gam,rate,AR);
[mesh.p,mesh.t, mesh.dgnodes] = mkmesh_ring(pde.porder,80,20,R0,R1,10);
% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1,:).^2+p(2,:).^2-R0^2)<1e-6, @(p) abs(p(1,:).^2+p(2,:).^2-R1^2)<1e-6};
mesh.boundarycondition = [1 2];  % Inner, Outer
% expressions for curved boundaries
mesh.curvedboundary = [1 1];
mesh.curvedboundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)-R0, @(p) sqrt(p(1,:).^2+p(2,:).^2)-R1};