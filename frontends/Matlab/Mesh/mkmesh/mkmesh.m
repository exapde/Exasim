function mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype,fb,fd,varargin)
%MKMESH Creates mesh data structure 
%   MESH=MKMESH(p,t,porder,bndexpr,elemtype,nodetype,fd,fdparams)
%
%      MESH:      Mesh structure
%      P:         Node positions
%      T:         Node connectivities
%      PORDER:    Polynomial Order of Approximation (default=1)
%      BNDEXPR:   Cell Array of boundary expressions. The 
%                 number of elements in BNDEXPR determines 
%                 the number of different boundaries
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution
%      FB:        Index for the curved boundary
%      FD:        Distance Function d(x,y,z)
%      FDPARAMS:  Additional parameters passed to FD
%
%   See also: MKT2F, SETBNDNBRS, MASTERNODES, CREATENODES
%

if nargin<4
    error('Require at least four input aguments');
end
if nargin<5, elemtype=0; end
if nargin<6, nodetype=0; end
if nargin<7, fb=[]; end
if nargin<8, fd=[]; end

%[p,t]=fixmesh(p,t);
dim = size(p,2);

mesh.porder = porder; % polynomial degree
mesh.p = p;  % node positions
mesh.t = t;  % element-to-node connectivities 

% compute face-to-node connectivities and element-to-face connectivities 
[mesh.f,mesh.t2f] = mkt2f(mesh.t,elemtype);

% make corrections for faces on the boundaries
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

% compute flags for curved faces and elements
if isempty(fb)
    mesh.fcurved = zeros(size(mesh.f,1),1); 
    mesh.tcurved = zeros(size(mesh.t,1),1);
else
    mesh.fcurved = (mesh.f(:,end)==-fb);
    ic = mesh.fcurved;
    mesh.tcurved = false(size(mesh.t,1),1);
    mesh.tcurved(mesh.f(ic,end-1)) = true;
end

% node positions on the master element
[mesh.plocal,mesh.tlocal,mesh.plocfc,mesh.tlocfc,permnode,permedge,permface] = mkmasternodes(mesh.porder,dim,elemtype,nodetype);
mesh.permnode = permnode; % positions of nodes at the vertices
mesh.permedge = permedge; % positions of nodes on the edges
mesh.permface = permface; % positions of nodes on the faces
% positions of nodes on the faces
if dim==1
    mesh.perm = permnode;
elseif dim==2
    mesh.perm = permedge;
elseif dim==3
    mesh.perm = permface;
end

% generate DG nodes
mesh.dgnodes = mkdgnodes(mesh,fd,varargin);

mesh.ne = size(mesh.t,1); % number of elements
mesh.np = size(mesh.p,1); % number of nodes 
mesh.nf = size(mesh.f,1); % number of faces
mesh.nd = dim;            % problem dimension
mesh.elemtype = elemtype; % 0 -> simplex elements or 1 -> tensor elements 
mesh.nodetype = nodetype; % 0 -> uniform distribution or 1 -> nonuniform distribution 
mesh.bndexpr = bndexpr;
mesh.fb = fb;
