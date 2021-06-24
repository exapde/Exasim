function mastergrid = mkmastersubgrid(pgeom,porder,pgauss,nref,nd,elemtype,nodetype)
%MASTERELEMENT  Create master grid structure
%    MASTER=MASTERGRID(PORDER,PGAUSS,ND,ELEMTYPE,NODETYPE,NREF)
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

% [AE,FE,DUDG,DUDG_DUH] = hdgassemble(UDG,UH,SH,mesh.dgnodes,mesh.bf,master.shapvt,master.shapvg,...
%         master.shapvgdotshapvl,master.shapft,master.shapfg,master.shapfgdotshapfc,...
%         master.perm-1,app.param,app.flags,app.factors,app.bcm,app.bcs',app.time);

mastergrid.nd     = nd;          % problem dimension
mastergrid.elemtype = elemtype;  % element type
mastergrid.nodetype = nodetype;  % node type
mastergrid.pgeom  = pgeom;       % geometry order
mastergrid.porder = porder;      % interpolation order 
mastergrid.pgauss = pgauss;      % integration order 
mastergrid.nref   = nref;        % refinement level 

%(pgeom,nd,elemtype,nodetype,pstar,nstar,pgstar,pdual,ndual,pgdual)

% standard master elements of order pgeom
mastergeom = masterelement(pgeom,2*pgeom,nd,elemtype,nodetype);

% standard master elements of order porder
masterelem = masterelement(porder,pgauss,nd,elemtype,nodetype);

% grid
[mastergrid.p,mastergrid.t] = masternodes(min(2*nref,1),nd,elemtype,nodetype);

% geometry nodes
mastergrid.geomnodes = mknodes(mastergrid.p,mastergrid.t,mastergeom.plocvl);

% interpolation nodes
mastergrid.elemnodes = mknodes(mastergrid.p,mastergrid.t,masterelem.plocvl);

% positions of geometry nodes on the faces
mastergrid.permgeom = mastergeom.perm;

% positions of interpolation nodes on the faces
mastergrid.permelem = masterelem.perm;

% volume shape functions at Gauss points on the standard master element
mastergrid.shapvl = masterelem.shapvl;

% volume shape functions at Gauss points times Gauss weights on the standard master element
for d=1:nd+1
    mastergrid.shapvg(:,:,d) = masterelem.shapvl(:,:,d)*diag(masterelem.gwvl);    
end

% face shape functions at Gauss points
mastergrid.shapfc = masterelem.shapfc;

% face shape functions at Gauss points times Gauss weights
for d=1:nd    
    mastergrid.shapfg(:,:,d) = masterelem.shapfc(:,:,d)*diag(masterelem.gwfc);
end

% shape functions of the standard master element at the geometry nodes
for j=1:size(mastergrid.geomnodes,3)
    mastergrid.shapgeom(:,:,:,j) = mkshape(pgeom,mastergeom.plocvl,mastergrid.geomnodes,elemtype);
end

% % Gauss points and weights on the master volume element  
% [master.gpvl,master.gwvl] = gaussquad(pgauss,nd,elemtype);
% 
% % shape functions and derivatives on the master volume element  
% master.shapvl = mkshape(porder,master.plocvl,master.gpvl,elemtype);
% 
% if nd>1
%     % Gauss points and weights on the master face element  
%     [master.gpfc,master.gwfc] = gaussquad(pgauss,nd-1,elemtype);
% 
%     % shape functions and derivatives on the master face element  
%     master.shapfc = mkshape(porder,master.plocfc,master.gpfc,elemtype);
% else
%     master.plocfc = 0;
%     master.gpfc   = 0;
%     master.gwfc   = 1;
%     master.shapfc = 1;
% end
% 
% % mass and convection matrices on the master element
% master.mass = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,1))';
% for ii=1:nd
%     master.conv(:,:,ii) = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,ii+1))';
% end   
% 
% master.ngv = size(master.gpvl,1);   % number of gasss points per element
% master.ngf = size(master.gpfc,1);   % number of gasss points per face
% master.npv = size(master.plocvl,1); % number of nodal points per element
% master.npf = size(master.plocfc,1); % number of nodal points per face
% 
% 
