function mesh=mkmesh_flatplate(ascale,nscale,gd,porder)

elemtype = 1;

% Parameters
if nargin<1, ascale=1; end
if nargin<2, nscale=1; end

% Geometry Parameters
H=0.1;   % chanel height
L=0.25;   % distance from freestream to leading edge of the plate
W=4.6870; % length of the plate

if nargin < 3, gd=1; end;
switch gd
    case 1    
    n1=7*nscale-1;
    n2=25*nscale-1;
    n3=19*nscale-1;     
    case 2
    n1=5*nscale-1;
    n2=21*nscale-1;
    n3=15*nscale-1;    
    case 3
    n1=2*nscale;
    n2=7*nscale-1;
    n3=5*nscale;    
    otherwise
    error('Unknown case.');
end

% Scaling Parameters
a1=2.2*ascale;
a2=3.2*ascale;
a3=5.1*ascale;


% Meshes in parameter space
[p01,t01]=squaremesh(n1,n3,1,elemtype); 
[p02,t02]=squaremesh(n2,n3,1,elemtype);

p01 = p01'; t01 = t01';
p02 = p02'; t02 = t02';

np01=size(p01,1); nt01=size(t01,1);
np02=size(p02,1); nt02=size(t02,1);

% Wall treatment
p01(:,1)=logdec(p01(:,1),a1);
p01(:,2)=loginc(p01(:,2),a3);
p02(:,1)=loginc(p02(:,1),a2);
p02(:,2)=loginc(p02(:,2),a3);

% Map
% p1=pmap(p01,[0,L,0,L],[0,0,H,H]); 
% p2=pmap(p02,[L,L+W,L,L+W],[0,0,H,H]);

p1=mapp(p01,[[0,L,0,L]' [0,0,H,H]']); 
p2=mapp(p02,[[L,L+W,L,L+W]' [0,0,H,H]']);

% Join
p=[p1;p2];
t=[t01;t02+np01];

% Remove duplicated nodes and orient elements
[p,t]=fixmesh(p,t);
p(:,1)=p(:,1)-L;

bndexpr =  {'all(p(:,2)<1e-6 & p(:,1)<1e-6)',...
            'all(p(:,2)<1e-6 & p(:,1)>-1e-6)',...
            'all(p(:,1)>4.6870-1e-6)', ...
            'all(p(:,2)>0.1-1e-6)',...
            'all(p(:,1)<-0.25+1e-6)'};     
        
mesh = mkmesh(p,t,porder,bndexpr,elemtype,1);        
mesh.p = mesh.p';
mesh.t = mesh.t';

xmin = min(mesh.p(1,:));
xmax = max(mesh.p(1,:));
ymin = min(mesh.p(2,:));
ymax = max(mesh.p(2,:));

mesh.boundaryexpr = {
  @(p) (p(1,:) < xmin + 1e-6), ... %  Freestream 
  @(p) (p(2,:) < ymin + 1e-6 && p(1,:) < 0.0 + 1e-6), ... %  Freestream  
  @(p) (p(2,:) < ymin + 1e-6 && p(1,:) > 0.0 - 1e-6), ... %  Wall  
  @(p) (p(1,:) > xmax - 1e-6),  ...%  Outflow
  @(p) (p(2,:) > ymax - 1e-6)         %  Freestream
};

mesh.periodicboundary = [];
mesh.periodicexpr = {};

mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.xpe = mesh.plocal; 
mesh.telem = mesh.tlocal;

figure(1);clf;boundaryplot(mesh,1); hold on; %Freestream
boundaryplot(mesh,2); %Freestream
boundaryplot(mesh,3); %wall
boundaryplot(mesh,4); %outflow
boundaryplot(mesh,5); %freestream 

mesh.xpe = mesh.plocal;
mesh.telem = mesh.tlocal;

figure(2); clf; meshplot(mesh);
