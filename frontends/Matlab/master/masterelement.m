function master = masterelement(porder,pgauss,nd,elemtype,nodetype)
%MASTERELEMENT  Create master element structure
%    MASTER=MASTERELEMENT(PORDER,PGAUSS,ND,ELEMTYPE,NODETYPE)
%
%      PORDER:    Polynomial Order of Approximation 
%      PGAUSS:    Degree of the polynomila to be integrated exactly
%                 (default: PGAUSS = 2*MESH.PORDER)
%      ND:        Dimension
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution
%

if nargin < 2, 
    pgauss = max(2*porder,1); 
end

master.nd     = nd;        % problem dimension
master.porder = porder;    % polynomial degree
master.pgauss = pgauss;    % order for Gauss points
master.elemtype = elemtype;% element type
master.nodetype = nodetype;% node type

% node positions on the master element
[master.plocvl,master.tlocvl,master.plocfc,master.tlocfc,master.permnode,...
    master.permedge,master.permface] = masternodes(porder,nd,elemtype,nodetype);

% positions of nodes on the faces
if nd==1
    master.perm = master.permnode;
elseif nd==2
    master.perm = master.permedge;
elseif nd==3
    master.perm = master.permface;
end

% Gauss points and weights on the master volume element  
[master.gpvl,master.gwvl] = gaussquad(pgauss,nd,elemtype);

% shape functions and derivatives on the master volume element  
master.shapvl = mkshape(porder,master.plocvl,master.gpvl,elemtype);

if nd>1
    % Gauss points and weights on the master face element  
    [master.gpfc,master.gwfc] = gaussquad(pgauss,nd-1,elemtype);

    % shape functions and derivatives on the master face element  
    master.shapfc = mkshape(porder,master.plocfc,master.gpfc,elemtype);
else
    master.plocfc = 0;
    master.gpfc   = 0;
    master.gwfc   = 1;
    master.shapfc = 1;
end

master.ngv = size(master.gpvl,1);   % number of gasss points per element
master.ngf = size(master.gpfc,1);   % number of gasss points per face
master.npv = size(master.plocvl,1); % number of nodal points per element
master.npf = size(master.plocfc,1); % number of nodal points per face


