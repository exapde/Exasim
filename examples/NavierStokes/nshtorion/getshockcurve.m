function xm = getshockcurve(mesh, sol)

[M, Mgrad] = gradient_machnumber(sol, 1.4);
x = mesh.dgnodes(:,1,:);
y = mesh.dgnodes(:,2,:);
ymin = min(y(:));
ymax = max(y(:));
xmax = max(x(:));

Mgradmax = 20;
ind = find(Mgrad>=Mgradmax);
x = x(ind);
y = y(ind);
ind = x<=-0.1 | y >= 2.8;
x = x(ind);
y = y(ind);

% figure(1); clf; scaplot(mesh, real(M),[]);
% hold on;
%plot(x,y,'o','Linewidth',2,'MarkerSize',6);

% P = polyfit(x(:), y(:), 12);
% s = linspace(min(x(:)), max(x(:)), 10000)';
% z = polyval(P, s);
% plot(s,z,'-','Linewidth',1);


P = polyfit(y(:), x(:), 8);
s = loginc(linspace(ymin, ymax, 1000)',2);
z = polyval(P, s);
ind = z <= xmax;
z = z(ind);
s = s(ind);

%plot(z,s,'-r','Linewidth',1);

xm = [z(:) s(:)];

