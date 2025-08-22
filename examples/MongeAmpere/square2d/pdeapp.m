
% solve Poisson equation with  homogeneous Neumann boundary condition
pdeapp_poisson; 

% get the velocity field as the gradient of the solution of the poisson equation
mesh.vdg = sol(:,2:3,:); 

% original mesh;
mesh.xpe = master.xpe;
mesh.telem = master.telem;
mesh.elemtype = master.elemtype;

% mesh1 is the adaptive mesh
mesh1 = mesh; 

% solve the transport equation with initial condition u(x,y,t=0) = x
mesh.udg = mesh.dgnodes(:,1,:);
pdeapp_transport;

% the x-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,1,:) = solt(:,:,:,end);

% solve the transport equation with initial condition u(x,y,t=0) = y
mesh.udg = mesh.dgnodes(:,2,:);
pdeapp_transport;

% the y-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,2,:) = solt(:,:,:,end);

% plot the velocity field
figure(1); clf; scaplot(mesh,sol(:,2,:),[],2); axis on; axis equal; axis tight;
figure(2); clf; scaplot(mesh,sol(:,3,:),[],2); axis on; axis equal; axis tight;

% plot the mesh density function
x = (mesh.dgnodes(:,1,:));
y = (mesh.dgnodes(:,2,:));
r = sqrt(x.^2+y.^2);
rho = 1 + a1*sech(a2*(r.^2-a^2));
figure(3); clf; scaplot(mesh,rho(:,1,:),[],1); axis on; axis equal; axis tight;

% plot the original mesh
figure(4); clf; meshplot(mesh,1); axis on; axis equal; axis tight;

% plot the adaptive mesh
figure(5); clf; meshplot(mesh1,1); axis on; axis equal; axis tight;

