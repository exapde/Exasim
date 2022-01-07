function checkudg3d(mesh,UDG,nz,opt,box)

if nargin<3
    opt = [];
end
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);

gam = 1.4;
ne = size(UDG,3);
porder = mesh.porder;

r = UDG(:,1,:);
rmin = min(r(:))
ind = find(r(:)==rmin);
[~,i2,~,i4] = ind2sub([(porder+1)^2 (porder+1) ne/nz nz],ind);
r = reshape(r,[(porder+1)^2 (porder+1) ne/nz nz]);
r = reshape(r(:,i2,:,i4),[(porder+1)^2 1 ne/nz]);
[i2 i4]
pause

ru = UDG(:,2,:);
ru = reshape(ru,[(porder+1)^2 (porder+1) ne/nz nz]);
ru = reshape(ru(:,i2,:,i4),[(porder+1)^2 1 ne/nz]);
rv = UDG(:,3,:);
rv = reshape(rv,[(porder+1)^2 (porder+1) ne/nz nz]);
rv = reshape(rv(:,i2,:,i4),[(porder+1)^2 1 ne/nz]);
rw = UDG(:,4,:);
rw = reshape(rw,[(porder+1)^2 (porder+1) ne/nz nz]);
rw = reshape(rw(:,i2,:,i4),[(porder+1)^2 1 ne/nz]);
rE = UDG(:,5,:);
rE = reshape(rE,[(porder+1)^2 (porder+1) ne/nz nz]);
rE = reshape(rE(:,i2,:,i4),[(porder+1)^2 1 ne/nz]);
u = ru./r;
v = rw./r;
w = rw./r;

figure(1); clf; scaplot(mesh,r,[0.4 1.2],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(r(:));
[~,i2] = max(r(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Density','Fontsize',16);
axis(box);

figure(2); clf; scaplot(mesh,rE,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(rE(:));
[~,i2] = max(rE(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Total energy','Fontsize',16);
axis(box);

pres = (gam-1)*(rE - 0.5*(ru.*u + rv.*v+ rw.*w));
u2 = sqrt(u.^2+v.^2+v.^2+w.^2);
mach = u2./sqrt(gam*pres./r);

figure(3); clf; scaplot(mesh,pres,[0.5 1.8],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(pres(:));
[~,i2] = max(pres(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Pressure','Fontsize',16);
axis(box);

figure(4); clf; scaplot(mesh,mach,[0 1.3],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(mach(:));
[~,i2] = max(mach(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Mach number','Fontsize',16);
axis(box);

figure(5); clf; scaplot(mesh,u,[-0.3 1.7],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Horizontal velocity','Fontsize',16);
axis(box);

