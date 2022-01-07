function mesh = mkmesh_naca(siz,porder,xd,thick)
%MKMESH_NACA Creates 2D mesh data structure for a NACA foil.
%   MESH=MKMESH_NACA(SIZ,PORDER,XD,THICK)
%
%      MESH:      Mesh structure
%      SIZ:       Desired nominal element size (defualt=0.04)
%      PORDER:    Polynomial Order of Approximation (default=1)
%      XD(4):     [x_min, x_max, y_min, y_max] for far field boundary
%                 (default=[-2,4,-2,2])
%      THICK:     Airfoil thickness in percentage (default=12)
%
%   See also: DISTMESH2D, FIXMESH (part of DISTMESH), MKT2F, SETBNDNBRS, 
%             UNIFORMLOCALPNTS, CREATENODES, FDNACA, NACA, FHSET
%
if nargin>0 & siz==0.0; siz=0.08; end
if nargin<1, siz=0.08; porder=1; end
if nargin<2 | isempty(xd), xd=[-2,4,-2,2]; end     % [x_min, x_max, y_min, y_max]
if nargin<3, thick=12; end     

% bounding box
box=1.1.*[xd(1),xd(3);xd(2),xd(4)];    
% fixed points at trailing edge
%xtrail=[0.8;0.85;0.9;0.95]; ytrail=naca(xtrail,thick);
xtrail=[0.8:siz:0.98]'; ytrail=naca(xtrail,thick);
% grid fixed points
fix=[xd(1),xd(3);xd(1),xd(4);xd(2),xd(3);xd(2),xd(4);... % domain corners
     0,0; 1,0;...                           % leading and trailing edges
     xtrail,ytrail; xtrail,-ytrail];        % extra points at trailing edge
[mesh.p,mesh.t] = distmesh2d(@fdnaca,@fhset,siz,box,fix,xd,thick);

[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr = {'all(sqrt((p(:,1)-0.5).^2+p(:,2).^2)<0.6)', ...
           'all(sqrt((p(:,1)-0.5).^2+p(:,2).^2)>0.6)'};     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

mesh.fcurved = (mesh.f(:,4)<0);
ic = find(mesh.fcurved);
mesh.tcurved = repmat(false,size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

mesh.porder = porder;
% [mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
% mesh.dgnodes = createnodes(mesh,@fdnaca,xd,thick);
[mesh.plocal,mesh.tlocal,mesh.plocfc,mesh.tlocfc,permnode,permedge,permface] = masternodes(mesh.porder,2,0,0);
mesh.permnode = permnode;
mesh.permedge = permedge;
mesh.permface = permface;
mesh.perm = permedge;

mesh.dgnodes = createnodes(mesh,@fdnaca,xd,thick);
mesh.ne = size(mesh.t,1);
mesh.np = size(mesh.p,1);
mesh.nf = size(mesh.f,1);
mesh.nd = 2;
mesh.elemtype = 0;
mesh.nodetype = 0;


function d=fdnaca(p,xd,thick)
%FDNACA Distance function for NACA foil inside a rectangular domain. 
%  D=FDNACA(P,XD,THICK)
%
%      P:         Node positions (Nx2)
%      XD(4):     [x_min, x_max, y_min, y_max] for far field boundary
%      THICK:     Airfoil thickness in percentage
%
%   See also:  NACA, DPOLY, DRECTANGLE, DDIFF
%
th=pi:-pi/200:pi/2;
xt = cos(th)+1; xt(end)=[];  yt=naca(xt,thick);  
xb=fliplr(xt); yb=-naca(xb,thick);
pv=[xt',yt';1,0;xb',yb'];
dfoil=dpoly(p,pv);                          % distance to foil

drec=drectangle(p,xd(1),xd(2),xd(3),xd(4));     % distance to domain edge

d=ddiff(drec,dfoil);

%--------------------------------------------------------------------------

function y=naca(x,thick)
%NACA Equation for a NACA profile. 
%  Y=NACA(X,THICK)
%
%      Y:         y-coordinate
%      X:         x-coordinate
%      THICK:     Airfoil thickness in percentage
%
y=5*0.01*thick*(0.29690*sqrt(x)-0.12600*x-0.35160*x.^2+0.28430*x.^3-0.10150*x.^4);

%--------------------------------------------------------------------------

function h=fhset(p,varargin)
%FHSET Mesh size function for NACA foil. 
%  H=FSET(P)
%
%      P:         Node positions (Nx2)
%      H:         Desired spacing
%
%   See also:  DCIRCLE
%

% --- Positions of the source points ---
s1x=0; s1y=0;       % leading edge
s2x=1; s2y=0;       % trailing edge


% --- Growth parameters ---
h0    = 0.001;       % initial edge lenght
ratio  = 0.02;       % growth ratio
radius = 0.05;         % 0 for a point source

% --- Distance functions ---
h1 = h0 + max(ratio*dcircle(p,s1x,s1y,radius),0);
h2 = h0 + max(ratio*dcircle(p,s2x,s2y,radius),0);
h12 = min(h1,h2);


h4 = h0 + max(ratio*dcircle(p,0.3,0,radius),0);
h5 = h0 + max(ratio*dcircle(p,0.6,0,radius),0);
h45 = min(h4,h5);

h=min(h12,h45);


