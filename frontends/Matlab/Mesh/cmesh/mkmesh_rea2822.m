gridNum=1;
porder=4;

% --- Mesh
% if gridNum==1
%    n1=18*porder+1; n2=20*porder+1; n3=16*porder+1; n4=4*porder+1;
% elseif gridNum==2
%    n1=18*porder+1; n2=16*porder+1; n3=22*porder+1; n4=4*porder+1;
% end
% [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
%                     [40, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1], ...
%                     [20, 10, 10, 10, 10, 10, 10]*10);

if gridNum==1
   n1=10*porder+1; n2=22*porder+1; n3=16*porder+1;
end
[x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                    [20, 5, 5, 5, 5, 10, 1, 1, 1, 1, 10], ...
                    [1, 1, 1, 1, 1, 1, 1]*4e5);

load rae2822.mat;
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




%save msh3dg.mat msh porder;
