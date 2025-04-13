% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% initialize pde structure and mesh structure
[pde,mesh] = initializeexasim();

% Define a PDE model: governing equations, initial solutions, and boundary conditions
pde.model = "ModelD";          % ModelC, ModelD, ModelW
pde.modelfile = "pdemodel";    % name of a file defining the PDE model

% Choose computing platform and set number of processors
pde.platform = "cpu";         % choose this option if NVIDIA GPUs are available
pde.mpiprocs = 1;             % number of MPI processors
pde.hybrid = 1;               % 0 -> LDG, 1 -> HDG
pde.debugmode = 0;
pde.nd = 1;

% Set discretization parameters, physical parameters, and solver parameters
pde.porder = 3;             % polynomial degree
pde.pgauss = 2*pde.porder;
pde.torder = 1;          % time-stepping order of accuracy
pde.nstage = 1;          % time-stepping number of stages
pde.tau = 1.0;              % DG stabilization parameter

T0 = 2000;     % temperature
kB = 8.617e-5; % Boltzmann Constant
R = 8.314;     % ideal gas constant
L = 2e-5;      % physical length 10 micrometer

% Diffusion coefficient
D0 = 5.64*1e-8;
H = 1.09;
D = D0*exp(-H/(kB*T0));

% reaction rate
% C0 = 7.45e4;
% k0  = 1.5;
% E = 1.95;
% K = k0*C0*exp(-E/(kB*T0));
K = 18;

tc = L*L/D;         % physical time 
Dstar = tc*D/(L*L); % dimensionless diffusion coefficient should be equal to 1
Kstar = K*tc;       % dimensionless reaction rate

Cb = 1e-3;          % Dirichlet boundary condition
pde.physicsparam = [Dstar Kstar Cb];       % unit thermal conductivity
pde.ppdegree = 1;          % degree of polynomial preconditioner
pde.RBdim = 0;
pde.tau = 1.0;           % DG stabilization parameter

tfstar = 1; % dimensionless final time (always 1)
dt = 1e-2;  % timestep size
nsteps = tfstar/dt;
pde.dt = dt*ones(1,nsteps);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = dt; % visualization timestep size
pde.linearsolvertol = 1e-12; % GMRES tolerance

% create a 1D grid from 0 to 1
[mesh.p,mesh.t] = linemesh(200);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1))<1e-8, @(p) abs(p(1)-1<1e-8)};

mesh.boundarycondition = [1;2]; %set boundary condition for each boundary

pde.gencode = 1;
% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% u = sol(:,1,:,end); % numerical solution
% 
% % steady-state exact solution
% x = linspace(0,1,1000);
% l = sqrt(Kstar/Dstar);
% uex = Cb*(exp(l*x)+exp(2*l-l*x))/(1+exp(2*l));
% 
% figure(1); clf; plot(mesh.dgnodes(:),u(:),x,uex);

return;

L0=20;

figure(1); clf; 
plot(L0*mesh.dgnodes(:),u1(:),'-','LineWidth',2); 
hold on;
plot(L0*mesh.dgnodes(:),u2(:),'-','LineWidth',2); 
plot(L0*mesh.dgnodes(:),u3(:),'-','LineWidth',2); 
plot(L0*mesh.dgnodes(:),u4(:),'-','LineWidth',2); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$\mu \mathrm{m}$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$\mathcal{C} (\mathrm{mol}/\mathrm{m}^3)$", 'interpreter', 'latex', 'FontSize', 26);
set(gca, 'XTick', [0:1:10]*2);
set(gca, 'YTick', [0:1:10]/1e4);
set(gca,'TickLength',[0.02, 0.02])
axis([0 L0 0 Cb])
box on;
grid on;
leg = legend({'$T = 1000 \mathrm{K}$','$T = 1200 \mathrm{K}$','$T = 1600 \mathrm{K}$','$T = 2000 \mathrm{K}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [40,10];
print -dpng hfoxidation1d.png


% for i = 1:100
%   u = sol(:,1,:,i);
%   figure(3); clf; plot(mesh.dgnodes(:),u(:))
% end
% v0 = (1-tanh(1000*(mesh.dgnodes-0.1)));
% v = 0*sol(:,1,:,:);
% k1 = Kstar/Dstar*k0;
% for i = 1:100
%   if i== 1
%     v(:,:,:,i) = v0 + dt*k1*sol(:,1,:,i)./(1 + exp(100-1e7*sol(:,1,:,i)));
%   else
%     v(:,:,:,i) = v(:,:,:,i-1) + dt*k1*sol(:,1,:,i)./(1 + exp(100-1e7*sol(:,1,:,i)));
%   end
% end
% u = v(:,:,:,end);
% figure(2); clf; plot(mesh.dgnodes(:),u(:));


