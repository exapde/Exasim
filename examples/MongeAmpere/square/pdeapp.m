% Add Exasim to Matlab search path
cdir = pwd(); ii = strfind(cdir, "Exasim");
run(cdir(1:(ii+5)) + "/install/setpath.m");

[mesh, rho, drhodx, drhody] = mkmesh_square(50, 4, 1);

mesh1 = radaptivity(mesh, rho, drhodx, drhody, 1e-4);




















% pde.nd = size(mesh.dgnodes, 2);
% pde.elemtype = mesh.elemtype;
% pde.porder = mesh.porder;      
% pde.nodetype = 1;
% pde.pgauss = 2*pde.porder;
% master = Master(pde);
% mesh.xpe = master.xpe;
% mesh.telem = master.telem;
% 
% x = (mesh.dgnodes(:,1,:));
% y = (mesh.dgnodes(:,2,:));
% r = sqrt(x.^2+y.^2);
% a1 = 10;
% a2 = 100;
% a = 0.25;
% rho = 1 + a1*sech(a2*(r.^2-a^2));
% drhodx = -(2*a1*a2*x.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;
% drhody = -(2*a1*a2*y.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;
% 
% rho1 = a1*sech(a2*(r.^2-a^2));
% drho1dx = -(2*a1*a2*x.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;
% drho1dy = -(2*a1*a2*y.*sinh(a2*(- a^2 + x.^2 + y.^2)))./cosh(a2*(- a^2 + x.^2 + y.^2)).^2;
% 
% rho2 = 2*(1+x).*(1+y);
% drho2dx = 2*(1+y);
% drho2dy = 2*(1+x);
% 
% rho = 1 + rho1.*rho2;
% drhodx = drho1dx.*rho2 + rho1.*drho2dx;
% drhody = drho1dy.*rho2 + rho1.*drho2dy;
% 
% % rho = 1 + rho1 + rho2;
% % drhodx = drho1dx + drho2dx;
% % drhody = drho1dy + drho2dy;
% 
% mesh1 = radaptivity(mesh, rho, drhodx, drhody, 1e-4);
% 
% 
% % rho1 = rho-1 + (1.5+0.5*x).*(1.5+0.5*y);
% % drho1dx = drhodx + 0.5*(1.5+0.5*y);
% % drho1dy = drhody + 0.5*(1.5+0.5*x);
% 
% 
% mesh1 = radaptivity(mesh, rho, drhodx, drhody, 1e-4);
% 
% mesh2 = radaptivity(mesh, rho1, drho1dx, drho1dy, 1e-4);
% 
% [~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
% jac = reshape(jac,[],1,size(mesh.dgnodes,3));
% solhm = meshdensity(mesh, sqrt(jac), 1e-1, 3);
% r = solhm(:,1,:);
% a = min(r(:));
% b = max(r(:));
% solhm = solhm/a;
% 
% rho1 = rho.* solhm(:,1,:); 
% drho1dx = drhodx.* solhm(:,1,:) - rho.* solhm(:,2,:); 
% drho1dy = drhody.* solhm(:,1,:) - rho.* solhm(:,3,:); 
% 
% r = rho1(:,1,:);
% a = min(r(:));
% b = max(r(:));
% rho1 = 1 + 5*(rho1-a)/(b-a);
% drho1dx = 5*(drho1dx-a)/(b-a);
% drho1dy = 5*(drho1dy-a)/(b-a);
% 
% mesh2 = radaptivity(mesh, rho1, drho1dx, drho1dy, 1e-4);
% 
% 
% % mesh0 = mesh;
% % mesh0.p(1,:) = inverselog(0.5*(1+mesh.p(1,:)),2);
% % mesh0.p(2,:) = inverselog(0.5*(1+mesh.p(2,:)),2);
% % mesh0.dgnodes(:,1,:) = inverselog(0.5*(1+mesh.dgnodes(:,1,:)),2);
% % mesh0.dgnodes(:,2,:) = inverselog(0.5*(1+mesh.dgnodes(:,2,:)),2);
% % figure(1); clf; meshplot(mesh0);
% % 
% % figure(2); clf; scaplot(mesh0, rho, [], 2, 1);
% % 
% % [~, ~, jac] = volgeom(master.shapent,permute(mesh.dgnodes,[1 3 2]));
% % figure(7); clf; scaplot(mesh, jac,[] , 2);
% % 
% % x = mesh.dgnodes(:,1,:);
% % y = mesh.dgnodes(:,2,:);
% % 
% % % b = 0;
% % % c = 2;
% % % rho1 = b + (rho - 1) .* exp(c*x) .* exp(c*y);
% % % drho1dx = drhodx .* exp(c*x) .* exp(c*y) + c*(rho - 1) .* exp(c*x) .* exp(c*y);
% % % drho1dy = drhody .* exp(c*x) .* exp(c*y) + c*(rho - 1) .* exp(c*x) .* exp(c*y);
% % 
% % rhomax = max(rho(:));
% % rhomin = min(rho(:));
% % rho1 = 0.5 + (rho - rhomin)/(rhomax - rhomin);
% % drho1dx = drhodx/(rhomax-rhomin);
% % drho1dy = drhody/(rhomax-rhomin);
% % 
% % mesh1 = radaptivity(mesh, rho, drhodx, drhody, 1e-4);
% % 
% % figure(6); clf; scaplot(mesh1, mesh1.vdg(:,1,:),[],2);
% % 
