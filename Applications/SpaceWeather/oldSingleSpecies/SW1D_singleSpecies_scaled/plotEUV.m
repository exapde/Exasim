% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/Installation/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
% pde.platform = "gpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;              % number of MPI processors

[pde,mesh] = pdeparams(pde,mesh);

%timeStep (s)
tstep = 15*60;
%Length of simulation (days)
tSimulation = 1;
%frequency of data (minutes)
freq = 15;

t0 = pde.physicsparam(19);
periodDay = 86400;
tstepStar = tstep/t0;
nTimeSteps = ceil(tSimulation*periodDay/tstep);
freqTimeSteps = ceil(freq*60/tstep);

pde.dt = tstepStar*ones(nTimeSteps,1);   % time step sizes
pde.visdt = pde.dt(1);         % visualization timestep size
pde.saveSolFreq = freqTimeSteps;          % solution is saved every 100 time steps
pde.soltime = freqTimeSteps:freqTimeSteps:length(pde.dt); % steps at which solution are collected

% search compilers and set options
pde = setcompilers(pde);       

% generate input files and store them in datain folder
[pde,mesh,master,dmd] = preprocessing(pde,mesh);

dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,pde.porder);    
nt = length(pde.soltime);
nElem = size(dgnodes,3);
nPoints = size(dgnodes,1);
sol = zeros([nPoints,1,nElem,nt]);
mesh.dgnodes = dgnodes;

UDG = initu(mesh,{0,0,0,0,0,0});
UDG = initialCondition1D(mesh,UDG,pde.physicsparam);
% load('IC0');

t = 0;
for iT = 1:nt
    indexT = pde.soltime(iT);
    t = t + pde.saveSolFreq*pde.dt(indexT);
    fprintf(('Solving for time %d\n'),t);
    for iElem = 1:nElem
        for iPoint = 1:nPoints
            xg = mesh.dgnodes(iPoint,:,iElem);
            u0 = UDG(iPoint,:,iElem);
            EUV0 = EUVsource1D(u0,xg,t,pde.physicsparam,pde.externalparam);
            sol(iPoint,:,iElem,iT) = [EUV0];
        end
    end
end

figure(1);clf;
for i=1:size(sol,4)
    EUV = sol(:,1,:,i); 
    
    plot(dgnodes(:),EUV(:),'linewidth',2);
    axis([min(dgnodes(:)) max(dgnodes(:)) 0 8e-3])
    grid on
    grid minor
    
%     text(min(dgnodes(:)) + 10,-11,sprintf('t = %d',i))

    pause(.3)
    hold off

end