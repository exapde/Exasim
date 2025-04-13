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

L = 4e-5;                % m physical length of the solid 
Twall = 1000;            % K wall temperature  
T0  = 1500;              % K  temperature used for nondimensionalization
R = 8.314;               % J/(K * mol) ideal gas constant
kB = 8.617e-5;           % eV/K Boltzmann Constant
k = 18;                  % 1/s reaction rate
H0 = 1.415*1e10;         % J/m^3  chemistry modulus
mu0 = -104.9 * 1e3;      % J/mol reference chemical potential 
C0 = 7.45 * 1e4;         % mol/m^3 reference concentration
D0 = 9.28*1e-6;          % m^2/s diffusion rate constant
H = 1.46;                % eV diffusion barrier energy
D1 = D0*exp(-H/(kB*T0)); % m^2/s diffusion rate with respect to T0
tc = L*L/D1;             % s physical final time  
Dc = tc*D0/(L*L);        % dimensionless diffusion constant

rhoHf = 13000;           % kg/m^3 Hf density
cpHf  = 170;             % J/(kg * K) Hf heat capacity
kappaHf = 15;            % J/(s * m * K) Hf thermal conductivity
CvHf = cpHf*rhoHf;       % J/(m^3 * K) Hf volumetric specific heat;
alphaHf = kappaHf/CvHf;  % m^2/s Hf thermal diffusivity
DTHf = tc*alphaHf/(L*L); % dimensionless thermal diffusivity

mu = [L Twall T0 R kB k H0 mu0 C0 D0 H tc Dc DTHf CvHf];

pde.physicsparam = mu;   % phyiscal parameters
pde.ppdegree = 1;        % degree of polynomial preconditioner
pde.RBdim = 0;
pde.tau = 1.0;           % DG stabilization parameter

tfstar = 1; % dimensionless final time (always 1)
dt = 1e-2/2;  % timestep size
nsteps = tfstar/dt;
pde.dt = dt*ones(1,nsteps);   % time step sizes
pde.soltime = 1:length(pde.dt); % steps at which solution are collected
pde.visdt = dt; % visualization timestep size
pde.linearsolvertol = 1e-12; % GMRES tolerance

% create a 1D grid from 0 to 1
[mesh.p,mesh.t] = linemesh(300);

% expressions for domain boundaries
mesh.boundaryexpr = {@(p) abs(p(1))<1e-8, @(p) abs(p(1)-1<1e-8)};

mesh.boundarycondition = [1;2]; %set boundary condition for each boundary

pde.gencode = 1;
% call exasim to generate and run C++ code to solve the PDE model
[sol,pde,mesh,master,dmd] = exasim(pde,mesh);

% for i = 1:length(pde.dt)
%   u = sol(:,1,:,i);
%   figure(3); clf; plot(mesh.dgnodes(:),u(:));
%   axis([0 1 0 1]);
%   pause(0.1);
% end

% u = sol(:,1,:,end); % numerical solution
% 
% % steady-state exact solution
% x = linspace(0,1,1000);
% l = sqrt(Kstar/Dstar);
% uex = Cb*(exp(l*x)+exp(2*l-l*x))/(1+exp(2*l));
% 
% figure(1); clf; 
% plot(mesh.dgnodes(:),u(:),'-',x,uex,'-','LineWidth',2);
% set(gca,'FontSize',22); 
% set(gca,'LooseInset',get(gca,'TightInset'))
% xlabel("$\mu \mathrm{m}$", 'interpreter', 'latex', 'FontSize', 26);
% ylabel("$\mathcal{C} (\mathrm{mol}/\mathrm{m}^3)$", 'interpreter', 'latex', 'FontSize', 26);
% % set(gca, 'XTick', [0:1:10]);
% % set(gca, 'YTick', [0:1:10]/1e4);
% % set(gca,'TickLength',[0.02, 0.02])
% % axis([0 10 0 Cb])
% % box on;
% % grid on;
% 
% % leg = legend({'$\mbox{GN}$','$\mbox{EIM-GN} (L=1, M=2N)$','$\mbox{EIM-GN} (L=1, M=3N)$','$\mbox{FOEIM-GN} (L=3, M=2N)$','$\mbox{FOEIM-GN} (L=3, M=3N)$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
% % leg.ItemTokenSize = [40,10];
% % print -dpng ex3fig1.png

return;

L0=40;

u1 = sol(:,1,:,5);
u2 = sol(:,1,:,20);
u3 = sol(:,1,:,60);
u4 = sol(:,1,:,120);

figure(1); clf; 
plot(L0*mesh.dgnodes(:),u1(:),'-','LineWidth',2); 
hold on;
plot(L0*mesh.dgnodes(:),u2(:),'-','LineWidth',2); 
plot(L0*mesh.dgnodes(:),u3(:),'-','LineWidth',2); 
plot(L0*mesh.dgnodes(:),u4(:),'-','LineWidth',2); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$\mu \mathrm{m}$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$\mathcal{E} = \mathcal{C}/\mathcal{C}_0$", 'interpreter', 'latex', 'FontSize', 26);
set(gca, 'XTick', [0:1:10]*6);
set(gca, 'YTick', [0:1:10]/10);
set(gca,'TickLength',[0.02, 0.02])
axis([0 L0 0 Cb])
box on;
grid on;
leg = legend({'$t = 0.004 \mbox{ } \mathrm{s}$','$t = 0.016 \mbox{ } \mathrm{s}$','$t = 0.048 \mbox{ } \mathrm{s}$','$t = 0.096 \mbox{ } \mathrm{s}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
%leg = legend({'$T = 1000 \mathrm{K}$','$T = 1200 \mathrm{K}$','$T = 1600 \mathrm{K}$','$T = 2000 \mathrm{K}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [40,10];
print -dpng hfoxidation1d.png




u1 = sol1000(:,1,:,81);
u2 = sol1500(:,1,:,80);
u3 = sol2000(:,1,:,120);

figure(1); clf; 
plot(L0*mesh.dgnodes(:),u1(:),'-','LineWidth',2); 
hold on;
plot(L0*mesh.dgnodes(:),u2(:),'-','LineWidth',2); 
plot(L0*mesh.dgnodes(:),u3(:),'-','LineWidth',2); 
set(gca,'FontSize',22); 
set(gca,'LooseInset',get(gca,'TightInset'))
xlabel("$\mu \mathrm{m}$", 'interpreter', 'latex', 'FontSize', 26);
ylabel("$\mathcal{E} = \mathcal{C}/\mathcal{C}_0$", 'interpreter', 'latex', 'FontSize', 26);
set(gca, 'XTick', [0:1:10]*6);
set(gca, 'YTick', [0:1:10]/10);
set(gca,'TickLength',[0.02, 0.02])
axis([0 L0 0 Cb])
box on;
grid on;
leg = legend({'$T = 1000 \mbox{ } \mathrm{K}, t = 0.160 \mbox{ } \mathrm{s}$','$T = 1500 \mbox{ } \mathrm{K}, t = 0.111 \mbox{ } \mathrm{s}$','$T = 2000 \mbox{ } \mathrm{K}, t = 0.048 \mbox{ } \mathrm{s}$','$t = 0.096 \mbox{ } \mathrm{s}$'}, 'interpreter', 'latex', 'FontSize', 20, 'Location', 'NE');
leg.ItemTokenSize = [40,10];
print -dpng hfoxidation1d_fig2.png


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


