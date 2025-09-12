function mesh = mkmesh_trefftz(m,n,porder,tparam)
%MKMESH_TREFFTZ Creates 2D mesh data structure for Trefftz airfoil.
%   MESH=MKMESH_TREFFTZ(M,N,PORDER,TPARAM)
%
%      MESH:      Mesh structure
%      M:         Number of points in the radial direction (default=15)
%      N:         Number of points in circumferential direction
%                 (default=30)
%      PORDER:    Polynomial Order of Approximation (default=3)
%      TPARAM:    Trefftz foil parameters
%                 TPARAM(1) = left x-shift of circle center 
%                             (trailing edge at (1,0)). (default=0.1)
%                 TPARAM(2) = y-shift of circle center. (default=0.05)
%                 TPARAM(3) = K-T exponent (=< 2) (2:Jukowski). (default=1.98)                      
%
%   See also: SQUAREMESH, FIXMESH (part of DISTMESH), MKT2F, SETBNDNBRS, 
%             UNIFORMLOCALPNTS, CREATENODES, TREFFTZ
%
% - Written by: J. Peraire
%

% Approx location of far field 
Rad = 20;
% Radial spacing (1,2, ... linear,quadratic ...)
rexp = 1.2;
elemtype=0;
nodetype=0;
dim=2;

if nargin<2, m=15; n=30; end
if nargin<3, porder=3; end
if nargin<4, tparam=[0.1,0.05,1.98]; end

n = 2*ceil(n/2);
[p0,t0] = squaremesh(m,n/2,0,elemtype);
[p1,t1] = squaremesh(m,n/2,1,elemtype);
np = size(p0,1);
t1 = t1+np;
p1(:,2) = p1(:,2)+1.0;
mesh.p = [p0; p1];
mesh.t = [t0; t1];
%[mesh.p,mesh.t]=fixmesh(mesh.p,mesh.t);
clear p0; clear p1; clear t0; clear t1;

mesh.porder = porder;
[mesh.plocal,mesh.tlocal,mesh.plocfc,mesh.tlocfc,permnode,permedge,permface] = masternodes(mesh.porder,dim,elemtype,nodetype);
mesh.permnode = permnode;
mesh.permedge = permedge;
mesh.permface = permface;
mesh.perm = permedge;

mesh.dgnodes = createnodes(mesh);

%First map to a cricle
mesh.p(:,1) = log(Rad)*mesh.p(:,1).^rexp;
mesh.p(:,2) = pi*mesh.p(:,2);

z = mesh.p(:,1) + i*mesh.p(:,2);
w = exp(z);
mesh.p = [real(w),imag(w)];
[mesh.p,mesh.t] = fixmesh(mesh.p, mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);
mesh.fcurved = repmat(true,size(mesh.f,1),1);
mesh.tcurved = repmat(true,size(mesh.t,1),1);

bndexpr = {'all(sqrt(sum(p.^2,2))<2)','all(sqrt(sum(p.^2,2))>2)'};  
% bndexpr = {'all(sqrt(sum(p.^2,2))<2)',...
%            'all(sqrt(sum(p.^2,2))>2) && any(p(:,1)<3+1e-6)',...
%            'all(sqrt(sum(p.^2,2))>2) && any(p(:,1)>=3-1e-6)'};  
%close all; simpplot(mesh.p,mesh.t);pause
% bndexpr = {'all(sqrt(sum(p.^2,2))<2)',...
%            'all(p(:,1)<5)',...
%            'all(p(:,1)>5)'};  
% bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-6 & p(:,1)<1/6+1e-6)',...
%            'all(p(:,2)<min(p0(:,2))+1e-6 & p(:,1)>=1/6-1e-6)',...
%            'all(p(:,1)>max(p0(:,1))-1e-6)', ...
%            'all(p(:,2)>max(p0(:,2))-1e-6)','all(p(:,1)<min(p0(:,1))+1e-6)'};     

mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

ind = find(abs(mesh.p(:,1).^2+mesh.p(:,2).^2-1)<1e-10);

%Now let's try a K-T transformation
x0 = tparam(1);
y0 = tparam(2);
n  = tparam(3);

%First Rotate to ensure that a point stays at the trailing edge
rot = atan2(y0,1+x0);
r = sqrt((1+x0)^2 + y0^2);

w = mesh.p(:,1) + i*mesh.p(:,2);
w = r*exp(-i*rot)*w + complex(-x0,y0);

%Now K-T
z = ((w-1)./(w+1)).^n;
w = ((1+z)./(1-z))*n;
mesh.p = [real(w),imag(w)];

% LE = mesh.p(ind,1)==min(mesh.p(ind,1)); LE = ind(LE);
% TE = mesh.p(ind,1)==max(mesh.p(ind,1)); TE = ind(TE);
% 
% chord = mesh.p(TE,1)-mesh.p(LE,1);
% mesh.p = mesh.p/chord;
% xle    = mesh.p(LE,1);
% mesh.p(:,1) = mesh.p(:,1)-xle;

%Now the same for the dgnodes
z = log(Rad)*mesh.dgnodes(:,1,:).^rexp + i*pi*mesh.dgnodes(:,2,:);
w = exp(z);

w = r*exp(-i*rot)*w + complex(-x0,y0);

z = ((w-1)./(w+1)).^n;
w = ((1+z)./(1-z))*n;

mesh.dgnodes(:,1,:) = real(w);
mesh.dgnodes(:,2,:) = imag(w);
% mesh.dgnodes = mesh.dgnodes/chord;
% mesh.dgnodes(:,1,:) = mesh.dgnodes(:,1,:)-xle;


mesh.ne = size(mesh.t,1);
mesh.np = size(mesh.p,1);
mesh.nf = size(mesh.f,1);
mesh.nd = dim;
mesh.elemtype = elemtype;
mesh.nodetype = nodetype;



