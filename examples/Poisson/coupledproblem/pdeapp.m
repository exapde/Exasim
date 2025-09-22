% clear exasim data from memory
clear pde mesh master dmd sol;

% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

% create pde and mesh for each PDE model
pdeapp1;
pdeapp2;

% figure(3); clf; 
% scaplot(mesh1,sol1(:,1,:),[],2,1); 
% hold on;
% scaplot(mesh2,sol2(:,1,:),[],2,1); 
% axis on; axis equal; axis tight;
% 
% pde{1}.modelnumber = 1;
% pde{2}.modelnumber = 2;
% % call exasim to generate and run C++ code to solve the PDE models
% sol = exasim(pde,mesh);
% 
% figure(4); clf; 
% scaplot(mesh1,sol{1}(:,1,:),[],2,1); 
% hold on;
% scaplot(mesh2,sol{2}(:,1,:),[],2,1); 
% axis on; axis equal; axis tight;

pde{1}.modelnumber = 1; % model number for pde1
pde{2}.modelnumber = 2; % model number for pde2

% set an interface condition for the interface of pde1
% the size of this array must be equal to the number of boundaries for pde1
% a nonzero element of this array corresponds to a particular fint of pdemodel1
mesh{1}.interfacecondition = [0;1;0;0]; 

% a list of fbou components that are used to impose interface condition
pde{1}.interfacefluxmap = [1]; 

% set an interface condition for the interface of pde2
% the size of this array must be equal to the number of boundaries for pde2
% a nonzero element of this array corresponds to a particular fint of pdemodel2
mesh{2}.interfacecondition = [0;0;0;1]; 

% a list of fbou components that are used to impose interface condition
pde{2}.interfacefluxmap = [1]; 

[sol,pde,mesh,master,dmd,compilerstr,runstr] = exasim(pde,mesh);

mesh{1}.porder = pde{1}.porder;
mesh{2}.porder = pde{2}.porder;
figure(5); clf; 
scaplot(mesh{1},sol{1}(:,1,:),[],2,1); 
hold on;
scaplot(mesh{2},sol{2}(:,1,:),[],2,1); 
axis on; axis equal; axis tight;


