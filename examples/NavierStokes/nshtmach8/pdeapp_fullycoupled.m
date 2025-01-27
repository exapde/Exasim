pdeapp1;
pdeapp2;

pde{1}.modelnumber = 1; % model number for pde1 (the navier-stokes equations)

% set boundary and interface conditions for pde1
% the size of this array must be equal to the number of boundaries for pde1
% each nonzero element of this array corresponds to a particular fbou of pdemodel1
mesh{1}.boundarycondition = [5;2;1]; 

% set an interface condition for the interface of pde1
% the size of this array must be equal to the number of boundaries for pde1
% a nonzero element of this array corresponds to a particular fint of pdemodel1
mesh{1}.interfacecondition = [1;0;0]; % only the first boundary has the interface condition

% a list of numerical flux components that are used to impose interface condition
pde{1}.interfacefluxmap = [4]; % this is 4 because the heat flux is used to impose the interface condition


pde{2}.modelnumber = 2; % model number for pde2 (the heat equation)

% set boundary and interface conditions for pde1
% the size of this array must be equal to the number of boundaries for pde2
% each nonzero element of this array corresponds to a particular fbou of pdemode2
mesh{2}.boundarycondition = [3;2;1];

% set an interface condition for the interface of pde2
% the size of this array must be equal to the number of boundaries for pde2
% a nonzero element of this array corresponds to a particular fint of pdemodel2
mesh{2}.interfacecondition = [1;0;0]; % only the first boundary has the interface condition

% a list of numerical flux components that are used to impose interface condition
pde{2}.interfacefluxmap = [1]; % this is 1 because the heat equation has only one numerical flux

[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

mesh1 = mesh{1};
mesh1.dgnodes = mesh1.dgnodes(:,:,dmd{1}{1}.elempart);
mesh2 = mesh{2};
mesh2.dgnodes = mesh2.dgnodes(:,:,dmd{2}{1}.elempart);
figure(4); clf; scaplot(mesh2, (Tref/Tinf)*sol{2}(:,1,:));
hold on;
scaplot(mesh1, (Tref/Tinf)*eulereval(sol{1}, 't',gam,Minf),[]);
set(gca,'FontSize',20); axis equal; axis tight;



