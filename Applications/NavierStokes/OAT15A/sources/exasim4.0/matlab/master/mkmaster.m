function master=mkmaster(mesh,pgauss)
%MKMASTER  Create master element structure
%    MASTER=MKMASTER(MESH)
%
%      MESH:      Mesh data structure
%      PGAUSS:    Degree of the polynomila to be integrated exactly
%                 (default: PGAUSS = 2*MESH.PORDER)
%

if nargin < 2 
    pgauss = max(2*mesh.porder,1); 
end

master.nd     = mesh.nd;       % problem dimension
master.porder = mesh.porder;   % polynomial degree for solution
if isfield(mesh,'morder')
    master.morder = mesh.morder;   % polynomial degree for geometry
else
    master.morder = mesh.porder;   % polynomial degree for geometry
end
% master.plocvl = mesh.plocal;   % nodal points on the master volume element
% master.plocfc = mesh.plocfc;   % nodal points on the master face element
% master.perm   = mesh.perm;     % locations of the nodal points on the faces of the master elements

dim      = mesh.nd;
elemtype = mesh.elemtype;      % element type flag: simplex or tensor
nodetype = mesh.nodetype;
porder   = mesh.porder;
morder   = master.morder;

% node positions on the master element
[master.plocvl,master.tlocvl,master.plocfc,master.tlocfc,permnode,permedge,permface] = masternodes(porder,dim,elemtype,nodetype);
if dim==1
    master.perm = permnode;
elseif dim==2
    master.perm = permedge;
elseif dim==3
    master.perm = permface;
end

% geometry nodes on the master element and face
[master.plocmv,master.tlocmv,master.plocmf,master.tlocmf,permnode,permedge,permface] = masternodes(morder,dim,elemtype,nodetype);
if dim==1
    master.permgeom = permnode;
elseif dim==2
    master.permgeom = permedge;
elseif dim==3
    master.permgeom = permface;
end

% Gauss points and weights on the master volume element  
[master.gpvl,master.gwvl] = gaussquad(pgauss,dim,elemtype);
% master.gpnvl = master.plocvl;
% master.gwnvl = (1/size(master.plocvl,1)) * ones(size(master.plocvl,1),1);

% shape functions and derivatives on the master volume element  
master.shapvl = mkshape(porder,master.plocvl,master.gpvl,elemtype);
% master.shapnvl = mkshape(porder,master.plocvl,master.plocvl,elemtype);

% shape functions and derivatives on the master volume element  
master.shapmv = mkshape(morder,master.plocmv,master.gpvl,elemtype);
master.shapmv = permute(master.shapmv,[2 1 3]);

if dim>1
    % Gauss points and weights on the master face element  
    [master.gpfc,master.gwfc] = gaussquad(pgauss,dim-1,elemtype);
%     master.gpnfc = master.plocfc;
%     master.gwnfc = (1/size(master.plocfc,1)) * ones(size(master.plocfc,1),1);

    % shape functions and derivatives on the master face element  
    master.shapfc = mkshape(porder,master.plocfc,master.gpfc,elemtype);
%     master.shapnfc = mkshape(porder,master.plocfc,master.plocfc,elemtype);

    % shape functions and derivatives on the master face element  
    master.shapmf = mkshape(morder,master.plocmf,master.gpfc,elemtype);
    master.shapmf = permute(master.shapmf,[2 1 3]);  
else
    master.plocfc = 0;
    master.gpfc   = 0;
    master.gwfc   = 1;
    master.shapfc = 1;
%     master.gpnfc   = 0;
%     master.gwnfc   = 1;
%     master.shapnfc = 1;
end

% mass and convection matrices on the master element
master.mass = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,1))';
for ii=1:dim
    master.conv(:,:,ii) = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,ii+1))';
end   

master.ngv = size(master.gpvl,1);   % number of gasss points per element
master.ngf = size(master.gpfc,1);   % number of gasss points per face
master.npv = size(master.plocvl,1); % number of nodal points per element
master.npf = size(master.plocfc,1); % number of nodal points per face

%--------------- update MASTER structure ----------------%
npv = master.npv;
ngv = master.ngv;
dim = master.nd;
master.shapvgdotshapvl  = zeros(npv*npv,ngv,dim+1);      
for d=1:dim+1 
    master.shapvt(:,:,d) = master.shapvl(:,:,d)';
    master.shapvg(:,:,d) = master.shapvl(:,:,d)*diag(master.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            master.shapvgdotshapvl((ii-1)*npv+jj,:,d) = master.shapvg(jj,:,d).*master.shapvl(ii,:,1);                    
        end
    end            
end

npf = master.npf;
ngf = master.ngf;
% face shape functions and their derivatives 
master.shapfgdotshapfc  = zeros(npf*npf,ngf,dim);      
for d=1:dim
    master.shapft(:,:,d) = master.shapfc(:,:,d)';
    master.shapfg(:,:,d) = master.shapfc(:,:,d)*diag(master.gwfc);
    for ii=1:npf
        for jj = 1:npf
            master.shapfgdotshapfc((ii-1)*npf+jj,:,d) = master.shapfg(jj,:,d).*master.shapfc(ii,:,1);                    
        end
    end            
end


