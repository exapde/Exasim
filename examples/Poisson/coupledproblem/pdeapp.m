% clear exasim data from memory
clear pde mesh master dmd sol;

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% create pde and mesh for each PDE model
pdeapp1;
pdeapp2;

figure(3); clf; 
scaplot(mesh1,sol1(:,1,:),[],2,1); 
hold on;
scaplot(mesh2,sol2(:,1,:),[],2,1); 
axis on; axis equal; axis tight;

pde{1}.modelnumber = 1;
pde{2}.modelnumber = 2;
% call exasim to generate and run C++ code to solve the PDE models
[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

mesh{1}.porder= pde{1}.porder;
mesh{2}.porder= pde{2}.porder;
figure(4); clf; 
scaplot(mesh{1},sol{1}(:,1,:),[],2,1); 
hold on;
scaplot(mesh{2},sol{2}(:,1,:),[],2,1); 
axis on; axis equal; axis tight;

pde{1}.modelnumber = 1; % model number for pde1
pde{2}.modelnumber = 2; % model number for pde2


% map the flux vector of pde2 to the flux vector of pde1
% the size of this array must be equal to ncu2
% each nonzero element of this array maps to a flux component of pde1
pde{1}.interfacefluxmap = [1]; 

% set an interface condition for each boundary of pde1
% the size of this array must be equal to the number of boundaries for pde1
% each nonzero element of this array corresponds to a particular fint of pdemodel1
mesh{1}.interfacecondition = [0;2;0;0]; 

% map the flux vector of pde1 to the flux vector of pde2
% the size of this array must be equal to ncu1
% each nonzero element of this array maps to a flux component of pde2
pde{2}.interfacefluxmap = [1]; 

% set an interface condition for each boundary of pde2
% the size of this array must be equal to the number of boundaries for pde2
% each nonzero element of this array corresponds to a particular fint of pdemodel2
mesh{2}.interfacecondition = [0;0;0;2]; 

[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

figure(5); clf; 
scaplot(mesh{1},sol{1}(:,1,:),[],2,1); 
hold on;
scaplot(mesh{2},sol{2}(:,1,:),[],2,1); 
axis on; axis equal; axis tight;



