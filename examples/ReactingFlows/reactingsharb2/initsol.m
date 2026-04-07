load('solp2.mat')
addpath(char(srcdir + "/Modeling/CNS5air/"));

R = 287;
gam = 1.4;
gam1 = gam-1;
Minf = pde.physicsparam(4);
T_ref  = 260.6;
rho_ref = 0.00235;
v_ref =  6921;
p_ref = R * rho_ref .* T_ref;

% !mppequil -T 260.6 -P 1.7576167e2 -s 2,3 air_5
% tm = init_reacting5_from_ideal(rho_ref, rho_ref*[v_ref; 0], p_ref, T_ref);

rho = sol(:,1,:);
v = sol(:,2:3,:)./sol(:,1,:);
ke = 0.5*(v(:,1,:).*v(:,1,:)+v(:,2,:).*v(:,2,:));
p = gam1*(sol(:,4,:) - rho.*ke);
T = gam*Minf^2 * p./rho;

[rho_phys, v_phys, T_phys, p_phys, e_phys, rhoE_phys] = computePhysicalStateFromNondim(rho, v, T, rho_ref, v_ref, T_ref);

figure(1); clf; scaplot(mesh, rho_phys,[],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
hold on; plot(mesh.dgnodes(:,1,1),mesh.dgnodes(:,2,1),'o');

figure(2); clf; scaplot(mesh, v_phys(:,1,:),[0 v_ref],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

figure(3); clf; scaplot(mesh, T_phys(:,1,:),[T_ref 15000],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
hold on; plot(mesh.dgnodes(:,1,1),mesh.dgnodes(:,2,1),'o');

figure(4); clf; scaplot(mesh, p_phys(:,1,:),[p_ref 1e5],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
hold on; plot(mesh.dgnodes(:,1,1),mesh.dgnodes(:,2,1),'o');

v1 = v_phys(:,1,:);
v2 = v_phys(:,2,:);
[npe, nc, ne] = size(sol);
rho_species = zeros(5,npe*ne);
rhoE = zeros(npe,1,ne);
for i = 1:npe*ne
  [i p_phys(i), T_phys(i), [v1(i) v2(i)]]
  info = equilibrate(p_phys(i), T_phys(i), [v1(i) v2(i)]);
  rho_species(:,i) = info.rho_species(:);
  rho(i) = info.rho;
  rhoE(i) = info.rhoE;
  % r1 = rho_phys(i);
  % p1 = p_phys(i);
  % T1 = T_phys(i);
  % u1 = [v1(i); v2(i)];  
  % [ta, tb, tc] = init_reacting5_from_ideal(r1, r1*u1, p1, T1);
  % max(abs(ta(:)-rho_species(:,i)))
  % max(abs(tb(:)-rhoE(i)))/max(abs(rhoE(i)))
  % max(abs(tc.Y(:)-info.Y(:)))
  % pause
  %[rho_species(:,i), rhoE(i), info] = init_reacting5_from_ideal(r1, r1*u1, p1, T1);
end
rho_species = permute(reshape(rho_species, [5 npe ne]), [2 1 3]);

figure(4); clf; scaplot(mesh, rho(:,1,:),[0 1e-2],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
% 
figure(5); clf; scaplot(mesh, sol(:,1,:)*rho_ref,[0 1e-2],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

figure(6); clf; scaplot(mesh, rhoE,[0 2.5e5],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

rhoE_phys = p_phys/gam1 + 0.5*(sol(:,1,:)*rho_ref).*(v_phys(:,1,:).*v_phys(:,1,:) + v_phys(:,2,:).*v_phys(:,2,:));
figure(7); clf; scaplot(mesh, rhoE_phys,[],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

figure(8); clf; scaplot(mesh, sol(:,4,:),[],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

figure(7); clf; scaplot(mesh, rho_species(:,5,:)./rho,[],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);

figure(7); clf; scaplot(mesh, T_phys,[],1);
colorbar; colormap('jet'); axis on; axis equal; axis tight; set(gca,'FontSize',16);
