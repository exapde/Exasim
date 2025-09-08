function [p, t] = cookmembranegrid(m,n,elemtype)

[p, t] = squaremesh(m,n,0,elemtype);
p = mapp(p',[[0,48,0,48]' [0,44,44,60]']);
p = p';




% [mesh.f,mesh.t2f] = mkt2f(mesh.t);
% 
% bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
%            'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<1e-3)'};     
% mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);
% 
% % map p
% 
% mesh.fcurved = zeros(size(mesh.f,1),1);
% mesh.tcurved = zeros(size(mesh.t,1),1);
% 
% mesh.porder = porder;
% [mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
% mesh.dgnodes = createnodes(mesh);
