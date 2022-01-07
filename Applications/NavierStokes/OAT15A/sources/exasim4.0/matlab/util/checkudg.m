function checkudg(mesh,UDG,opt,box)

if nargin<3
    opt = [];
end
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);

u = UDG(:,1,:);
figure(1); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Density','Fontsize',16);
axis(box);

u = UDG(:,4,:);
figure(2); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Total energy','Fontsize',16);
axis(box);

u = eulereval(UDG(:,1:4,:),'p',1.4);
figure(3); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Pressure','Fontsize',16);
axis(box);

u = eulereval(UDG(:,1:4,:),'t',1.4,1);
figure(4); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Temperature','Fontsize',16);
axis(box);

u = eulereval(UDG(:,1:4,:),'h',1.4);
figure(5); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Enthalpy','Fontsize',16);
axis(box);

u = eulereval(UDG(:,1:4,:),'u',1.4);
figure(6); clf; scaplot(mesh,u,[],1,opt); axis tight; axis off; colormap jet; 
[~,i1] = min(u(:));
[~,i2] = max(u(:));
ii = [i1 i2];
hold on;
plot(x(ii),y(ii),'ok','MarkerSize',10);
title('Velocity','Fontsize',16);
axis(box);

