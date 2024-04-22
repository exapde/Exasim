gridNum=1;
porder=4;

if gridNum==1
   n1=24+1; n2=18+1; n3=24+1;
end
[x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                    [5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 5], ...
                    [1, 1, 1, 1, 1, 1, 1]*1e2);

thick = 12;
th = (pi:-pi/200:pi/2)';
xt = (cos(th)+1)*1.0089304129;  
xt = xt(end:-1:1);
yt=naca(xt,thick);  
xb = flipud(xt);   
yb=-naca(xb,thick);
xf =[xt; xb(2:end)];
yf =[yt; yb(2:end)];
xf(end) = xf(1);
yf(end) = yf(1);

[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

% create DG mesh
bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<1','true'};
msh=cart2msh(porder,xm,ym,bndexpr);
   
mesh.p = msh.p';
mesh.t = msh.t'+1;
[mesh.p,mesh.t] = fixmesh(mesh.p, mesh.t);
mesh.porder = porder;
%[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.plocal = msh.s;
mesh.tlocal = double(msh.tlocal);
mesh.dgnodes = msh.p1;
[mesh.f,mesh.t2f] = mkt2f(mesh.t);
mesh.fcurved = repmat(true,size(mesh.f,1),1);
mesh.tcurved = repmat(true,size(mesh.t,1),1);

mesh.f = setbndnbrs(mesh.p,mesh.f,msh.bndexpr);

figure(1); clf;
meshplot(mesh,1);
hold on;
plot(xf,yf,'r*');

