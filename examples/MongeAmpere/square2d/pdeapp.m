
% solve Poisson equation with  homogeneous Neumann boundary condition
pdeapp_poisson; % return sol 

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
mesh1.dgnodes(:,1,:) = solt(:,1,:,end);

% solve the transport equation with initial condition u(x,y,t=0) = y
mesh.udg = mesh.dgnodes(:,2,:);
pdeapp_transport;

% the y-component of the adaptive mesh is the solution of the transport
% equation at time t = 1
mesh1.dgnodes(:,2,:) = solt(:,1,:,end);

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

d = mesh1.dgnodes - mesh.dgnodes;
% plot the adaptive mesh
figure(6); clf; scaplot(mesh,d(:,1,:),[],1); axis on; axis equal; axis tight;
figure(7); clf; scaplot(mesh,d(:,2,:),[],1); axis on; axis equal; axis tight;

v1 = sol(:,2,:);
v2 = sol(:,3,:);
x1 = mesh.dgnodes(:,1,:);
x2 = mesh.dgnodes(:,2,:);
r = sqrt(x1.^2 + x2.^2);
rho = 1 + a1*sech(a2*(r.^2 - a^2));
F = rho/theta;
t = 0;
w1 = v1./(t + (1-t)*F);
w2 = v2./(t + (1-t)*F);

figure(6); clf; scaplot(mesh,w1,[],1); axis on; axis equal; axis tight;
figure(7); clf; scaplot(mesh,w2,[],1); axis on; axis equal; axis tight;
figure(8); clf; scaplot(mesh,F,[],1); axis on; axis equal; axis tight;

% Compute rho, drhodx, and drhody
% 0. Calculate rho from the flow field 
% 0.a gradient of the flow density = sqrt(UDG(:,5,:).^2 + UDG(:,9,:).^2);
% 0.b rho = 1 + limiting(grad r, 0, 10, 1e3, 0)
% 0.c Calculate drhodx and drhody

% INPUT: mesh, rho, drhodx, drhody
% OUTPUT: mesh1 
% mesh1 = radaptivity(mesh, rho, drhodx, drhody, diffcoeff)

% Implementation of radaptivity
% 0. mesh1 = mesh
% 1. Calculate theta = int_\Omega rho dx / int_\Omega dx
% 2. calculate F = rho/theta 
% 3. mesh.vdg(:,1,:) = F, mesh.vdg(:,2,:) = drhodx/theta, mesh.vdg(:,3,:) = drhody/theta
% 4. Solve the Poisson equation by calling pdeapp_poisson to obtain sol = (u, dudx, dudy)
% 5. mesh.vdg(:,4,:) = sol(:,2,:), mesh.vdg(:,5,:) = sol(:,3,:)
% 6. Solve the transport equation with IC u = x 
% 6.a mesh.udg = mesh.dgnodes(:,1,:) 
% 6.b call pdeapp_transport and set mesh1.dgnodes(:,1,:) = solt(:,1,:,end)
% 7. Solve the transport equation with IC u = y 
% 7.a mesh.udg = mesh.dgnodes(:,2,:) 
% 7.b call pdeapp_transport and set mesh1.dgnodes(:,2,:) = solt(:,2,:,end)

