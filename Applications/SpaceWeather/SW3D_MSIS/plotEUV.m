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
pde.mpiprocs = 4;              % number of MPI processors
polarCoordinates = 1;

[pde,mesh] = pdeparams(pde,mesh);

%timeStep (s)
tstep = 24*4*15*60;
%Length of simulation (days)
tSimulation = 0.25;
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

nt = length(pde.soltime);
nElem = size(mesh.dgnodes,3);
nPoints = size(mesh.dgnodes,1);
sol = zeros([nPoints,2,nElem,nt]);

% UDG = initu(mesh,{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
% UDG = initialCondition(mesh,UDG,pde.physicsparam);
% load('IC0');

UDG = mesh.udg;
r = UDG(:,1,:);
rho = exp(r);
T = UDG(:,5,:)./sqrt(rho);
sol(:,1,:) = log10(rho);
sol(:,2,:) = T;
% 
% t = 0;
% for iT = 1:nt
%     indexT = pde.soltime(iT);
%     t = t + pde.saveSolFreq*pde.dt(indexT);
%     fprintf(('Solving for time %d of %d\n'),iT,nt);
%     for iElem = 1:nElem
%         for iPoint = 1:nPoints
%             xg = mesh.dgnodes(iPoint,:,iElem);
%             u0 = UDG(iPoint,:,iElem);
%             EUV0 = EUVsource3D(u0,xg,t,pde.physicsparam,pde.externalparam);
%             sol(iPoint,3,iElem,iT) = EUV0;
%         end
%     end
% end

if polarCoordinates
    dgnodesPolar = zeros(size(mesh.dgnodes));
   [a,e,r] = cart2sph(mesh.p(1,:),mesh.p(2,:),mesh.p(3,:));
   p2 = [a' e' r']';
   
   ne = pde.ne;
   for ie = 1:ne
       dge = mesh.dgnodes(:,:,ie);
       [adg,edg,rdg] = cart2sph(dge(:,1),dge(:,2),dge(:,3));
       adg = mod(adg+pi,2*pi)-pi;
       dgnodesPolar(:,:,ie) = [adg edg rdg];
       
       te = mesh.t(:,ie);
       pe = p2(:,te);
       pe(1,:) = mod(pe(1,:)+pi,2*pi)-pi;

       sa = unique((pe(1,:)>0));
       if (length(sa)>1 && max(abs(pe(1,:))>2))
           p2(1,te) = pe(1,:) + (pe(1,:)<0)*2*pi;
           dgnodesPolar(:,1,ie) = adg + (adg<0)*2*pi;
       end
   end
   
   dgnodesPolar(:,1,:) = dgnodesPolar(:,1,:)*180/pi;
   dgnodesPolar(:,2,:) = dgnodesPolar(:,2,:)*180/pi;
   dgnodesPolar(:,3,:) = (dgnodesPolar(:,3,:)-pde.physicsparam(16))*pde.physicsparam(18)/1000;

   p2(1,:) = p2(1,:)*180/pi;
   p2(2,:) = p2(2,:)*180/pi;
   p2(3,:) = (p2(3,:)-pde.physicsparam(16))*pde.physicsparam(18)/1000;
   
   mesh.dgnodes = dgnodesPolar;
   mesh.p = p2;

end

pde.visscalars = {"density",1,"Temperature",2};  % list of scalar fields for visualization
pde.visvectors = {}; % list of vector fields for visualization
xdg = vis(sol,pde,mesh); % visualize the numerical solution