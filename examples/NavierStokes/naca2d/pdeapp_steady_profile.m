% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

porder = 4;                     % polynomial degree
tau = 3.0;                      % stabilization parameter
gam = 1.4;                      % gas constant
Minf = 0.2;                     % freestream mach number
alpha = 0*pi/180;               % angle of attack
rinf = 1.0;                     % freestream density
ruinf = cos(alpha);             % freestream horizontal velocity
rvinf = sin(alpha);             % freestream vertical velocity
pinf = 1/(gam*Minf^2);          % freestream pressure
rEinf = 0.5+pinf/(gam-1);       % freestream energy
Re = 100;                       % Reynolds number 
Pr = 0.72;                      % Prandtl number 
ui = [ 1, cos(alpha), sin(alpha), 0.5+pinf/(gam-1)];

% initialize pde structure and mesh structure
[pde,~] = initializeexasim();

pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.mpiprocs = 1;              % number of MPI processors
pde.hybrid = 1;
pde.debugmode = 0;
pde.porder = porder;
pde.pgauss = 2*porder;

pde.physicsparam = [gam Re Pr Minf rinf ruinf rvinf rEinf];
pde.tau = tau;              % DG stabilization parameter
pde.GMRESrestart = 400;
pde.linearsolvertol = 1e-7; % GMRES tolerance
pde.linearsolveriter = 400;
pde.RBdim = 0;
pde.neb = 512;

% naca mesh
mesh = mkmesh_naca0012(porder,1,1);

pde.platform = "gpu";
pde.GMRESortho = 1;
pde.ppdegree = 20;
[sol,pde,mesh] = exasim(pde,mesh);
return;

platforms  = ["cpu", "gpu"];
GMRESortho = [0,1];
ppdegree = 0:20;

for i = 1:2
   for j = 1:2
       for k = 1:21
           pde.platform = platforms(i);         % choose this option if NVIDIA GPUs are available
           pde.GMRESortho = GMRESortho(j);
           pde.ppdegree = ppdegree(k);


           % Create a file name using sprintf
           filename = sprintf('output_steady/platform_%s_GMRESortho_%d_ppdegree_%d.txt', pde.platform, pde.GMRESortho, pde.ppdegree);
           disp(filename)


           diary(filename)
           diary on


           % call exasim to generate and run C++ code to solve the PDE model
           [sol,pde,mesh] = exasim(pde,mesh);


       end
   end
end

% % plot solution
% mesh.porder = porder;
% mesh.dgnodes = createdgnodes(mesh.p,mesh.t,mesh.f,mesh.curvedboundary,mesh.curvedboundaryexpr,porder);    
% figure(1); clf; scaplot(mesh,eulereval(sol(:,1:4,:),'M',gam),[],2); axis off; axis equal; axis tight;
