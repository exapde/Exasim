function a = avf(mesh, master, app, UDG)


lambda = app.lambda;
S0 = app.S0;
kappa = app.kappa;    
alpha=1e3;
div = divergence(UDG,1);
% he = meshsize(mesh);
% hk =  (kappa*min(he(:)))^2;
divmax = max(div(:));
% hk = kappa^2*2e-4;
% hk = kappa^2*3e-3;
% hk = kappa^2*8e-4 * 2 /app.porder; %THISISTHEONETHATWORKS

hk = kappa^2*5e-5 *2/app.porder;
s = limiting(div,0,10,alpha,0);
% s = limiting(div,0,10,alpha,0);
% mesh.ib = [];
% mesh.in = 1:size(mesh.p2,1);
mesh.ib = [];
mesh.in = 1:size(mesh.p2,2);
s = cgpoisson(mesh, master, s, [hk 1.0]);        
s = s/max(s(:));
a = (s-S0).*(atan(alpha*(s-S0))/pi + 0.5) - atan(alpha)/pi + 0.5;    

a = reshape(a(mesh.t2), size(UDG,1), 1, mesh.ne);    
%[min(a(:)) max(a(:))]
dist = tanh(mesh.dist*5);
a = lambda*(a.*dist);         
%a = 0.01*dist;
if size(mesh.dgnodes,2)==3
mesh.dgnodes(:,3,:) = [];
end
% figure(1);clf;scaplot(mesh,a)
% figure(2);clf;scaplot(mesh,limiting(div,0,divmax/2,alpha,0))
% figure(3);clf;scaplot(mesh,reshape(s(mesh.t2'), size(UDG,1), 1, mesh.ne))
% pauseopen 