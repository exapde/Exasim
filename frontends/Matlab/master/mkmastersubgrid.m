function mastersubgrid = mkmastersubgrid(pgeom,porder,pgauss,nref,nd,elemtype,nodetype)
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

mastersubgrid.nd     = nd;          % problem dimension
mastersubgrid.elemtype = elemtype;  % element type
mastersubgrid.nodetype = nodetype;  % node type
mastersubgrid.pgeom  = pgeom;       % geometry order
mastersubgrid.porder = porder;      % interpolation order 
mastersubgrid.pgauss = pgauss;      % integration order 
mastersubgrid.nref   = nref;        % refinement level 

%(pgeom,nd,elemtype,nodetype,pstar,nstar,pgstar,pdual,ndual,pgdual)

% standard master elements of order pgeom
mastergeom = masterelement(pgeom,2*pgeom,nd,elemtype,nodetype);

% standard master elements of order porder
masterelem = masterelement(porder,pgauss,nd,elemtype,nodetype);

% standard master elements of order porder
mastergmlm = masterelement(porder,2*pgeom,nd,elemtype,nodetype);

% grid over the super-volume master element and super-face master element
[mastersubgrid.p,mastersubgrid.t,mastersubgrid.pf,mastersubgrid.tf] = masternodes(max(2^nref,1),nd,elemtype,nodetype);

% geometry nodes over the super-volume master element
mastersubgrid.geomnodes = mknodes(mastersubgrid.p,mastersubgrid.t,mastergeom.plocvl);

% interpolation nodes over the super-volume master element
mastersubgrid.elemnodes = mknodes(mastersubgrid.p,mastersubgrid.t,masterelem.plocvl);

% gauss nodes over the super-face master element
%mastersubgrid.gasfnodes = mknodes(mastersubgrid.pf,mastersubgrid.tf,masterelem.gpfc);
mastersubgrid.gasfnodes = mknodes(mastersubgrid.pf,mastersubgrid.tf,mastergmlm.gpfc);
    
% positions of geometry nodes on the faces
mastersubgrid.permgeom = mastergeom.perm;

% positions of interpolation nodes on the faces
mastersubgrid.permelem = masterelem.perm;
mastersubgrid.perm     = masterelem.perm;

% volume shape functions at Gauss points on the standard master element
mastersubgrid.shapvl = masterelem.shapvl;

% geometry shape functions at Gauss points on the standard master element
shapmv = mkshape(pgeom,mastergeom.plocvl,masterelem.gpvl,elemtype);
for i=1:size(shapmv,3)
    mastersubgrid.shapmv(:,:,i) = shapmv(:,:,i)';
end
    
% Gauss points and weights on the master volume element  
mastersubgrid.gpvl = masterelem.gpvl;
mastersubgrid.gwvl = masterelem.gwvl;

% node positions on the master volume element
mastersubgrid.plocvl = masterelem.plocvl;
mastersubgrid.tlocvl = masterelem.tlocvl;

% volume shape functions at Gauss points times Gauss weights on the standard master element
for d=1:nd+1
    mastersubgrid.shapvg(:,:,d) = masterelem.shapvl(:,:,d)*diag(masterelem.gwvl);    
end

% geometry shape functions at Gauss points on the master face element
shapmf = mkshape(pgeom,mastergeom.plocfc,masterelem.gpfc,elemtype);
for i=1:size(shapmf,3)
    mastersubgrid.shapmf(:,:,i) = shapmf(:,:,i)';
end

shapmh = mkshape(pgeom,mastergeom.plocfc,mastergmlm.gpfc,elemtype);
for i=1:size(shapmh,3)
    mastersubgrid.shapmh(:,:,i) = shapmh(:,:,i)';
end

% face shape functions at Gauss points
mastersubgrid.shapfc = masterelem.shapfc;

% Gauss points and weights on the master face element  
mastersubgrid.gpfc = masterelem.gpfc;
mastersubgrid.gwfc = masterelem.gwfc;

% node positions on the super master face element
mastersubgrid.plocfc = masterelem.plocfc;
mastersubgrid.tlocfc = masterelem.tlocfc;

% face shape functions at Gauss points times Gauss weights
for d=1:nd    
    mastersubgrid.shapfg(:,:,d) = masterelem.shapfc(:,:,d)*diag(masterelem.gwfc);
end

% geometry shape functions of the super-volume master element at the geometry nodes
for j=1:size(mastersubgrid.geomnodes,3)
    tmp = mkshape(pgeom,mastergeom.plocvl,mastersubgrid.geomnodes(:,:,j),elemtype);
    mastersubgrid.shapgeom(:,:,j) = tmp(:,:,1)';
end

% face shape functions of the super master element at the Gauss points  
for j=1:size(mastersubgrid.gasfnodes,3)
    tmp = mkshape(pgeom,mastergeom.plocfc,mastersubgrid.gasfnodes(:,:,j),elemtype);
    mastersubgrid.shapsf(:,:,j) = tmp(:,:,1)';
    mastersubgrid.shapsg(:,:,j) = mastersubgrid.shapsf(:,:,j)'*diag(mastergmlm.gwfc);
end

npv = size(mastersubgrid.shapvl,1);
ngv = size(mastersubgrid.shapvl,2);
mastersubgrid.shapvgdotshapvl  = zeros(npv*npv,ngv,nd+1);      
for d=1:nd+1
    mastersubgrid.shapvt(:,:,d) = mastersubgrid.shapvl(:,:,d)';
    mastersubgrid.shapvg(:,:,d) = mastersubgrid.shapvl(:,:,d)*diag(mastersubgrid.gwvl);    
    for ii=1:npv
        for jj = 1:npv
            mastersubgrid.shapvgdotshapvl((ii-1)*npv+jj,:,d) = mastersubgrid.shapvg(jj,:,d).*mastersubgrid.shapvl(ii,:,1);                    
        end
    end            
end

% face shape functions and their derivatives 
npf = size(mastersubgrid.shapfc,1);
ngf = size(mastersubgrid.shapfc,2);
nqf = size(mastergeom.shapfc,1);
nsf = size(mastersubgrid.tf,1);
nfe = size(mastersubgrid.perm,2);
npm = size(mastergeom.shapvl,1);

mastersubgrid.shapfgdotshapfc  = zeros(npf*npf,ngf,nd);   
mastersubgrid.shaphgdotshaphc  = zeros(npf*npf,size(mastergeom.shapfc,2),nd);   
for d=1:nd
    mastersubgrid.shapft(:,:,d) = mastersubgrid.shapfc(:,:,d)';
    mastersubgrid.shapfg(:,:,d) = mastersubgrid.shapfc(:,:,d)*diag(mastersubgrid.gwfc);
    mastersubgrid.shapht(:,:,d) = mastergmlm.shapfc(:,:,d)';
    mastersubgrid.shaphg(:,:,d) = mastergmlm.shapfc(:,:,d)*diag(mastergmlm.gwfc);
    for ii=1:npf
        for jj = 1:npf
            mastersubgrid.shapfgdotshapfc((ii-1)*npf+jj,:,d) = mastersubgrid.shapfg(jj,:,d).*mastersubgrid.shapfc(ii,:,1);                    
            mastersubgrid.shaphgdotshaphc((ii-1)*npf+jj,:,d) = mastersubgrid.shaphg(jj,:,d).*mastergmlm.shapfc(ii,:,1);                    
        end
    end            
end

for d=1:size(mastersubgrid.shapsf,3)
    for ii=1:size(mastersubgrid.shapsf,2)
        for jj = 1:npf            
            %mastersubgrid.shapfgdotshapsf((ii-1)*npf+jj,:,d) = mastersubgrid.shapfg(jj,:,1).*(mastersubgrid.shapsf(:,ii,d)');                    
            mastersubgrid.shaphgdotshapsf((ii-1)*npf+jj,:,d) = mastersubgrid.shaphg(jj,:,1).*(mastersubgrid.shapsf(:,ii,d)');                    
        end
    end            
end

for d=1:size(mastersubgrid.shapsf,3)
    for ii=1:size(mastersubgrid.shapsf,2)
        for jj = 1:size(mastersubgrid.shapsf,2)            
            %mastersubgrid.shapsgdotshapsf((ii-1)*nqf+jj,:,d) = (mastersubgrid.shapsf(:,jj,d)').*(mastersubgrid.shapsf(:,ii,d)')*diag(masterelem.gwfc);                    
            mastersubgrid.shapsgdotshapsf((ii-1)*nqf+jj,:,d) = (mastersubgrid.shapsf(:,jj,d)').*(mastersubgrid.shapsf(:,ii,d)')*diag(mastergmlm.gwfc);                    
        end
    end            
end

[mastersubgrid.f,mastersubgrid.t2f] = mkt2f(mastersubgrid.t,mastersubgrid.elemtype);
mastersubgrid.bf = reshape(mastersubgrid.f(abs(mastersubgrid.t2f'),end),[size(mastersubgrid.perm,2) size(mastersubgrid.t,1)]);
ind = find(mastersubgrid.bf==0);
mastersubgrid.bf(ind) = mastersubgrid.bf(ind) - 1;
mastersubgrid.f(:,end+1) = 0;
mastersubgrid.elcon = elconnectivities(mastersubgrid);
mastersubgrid.f(:,end) = [];

[p1,~,p1fc,~,permnode,permedge,permface] = masternodes(1,nd,elemtype,nodetype);
if nd==1
    perm1 = permnode;
elseif nd==2
    perm1 = permedge;
elseif nd==3
    perm1 = permface;
end

boudnodes = zeros(size(p1fc,1),nd,nsf,nfe);
for i=1:size(perm1,2)
    for j=1:size(mastersubgrid.tf,1)
        plocal = mastersubgrid.pf(mastersubgrid.tf(j,:),:);        
        ptm = getbounodes(p1(perm1(:,i),:),plocal);
        boudnodes(:,:,j,i) = getbounodes(ptm,p1fc);
    end
end

elemnodes = mknodes(mastersubgrid.p,mastersubgrid.t,p1);
mastersubgrid.bf(:,:,2:3) = 0;
for k=1:size(mastersubgrid.bf,2)
    for m=1:size(mastersubgrid.bf,1)
        if mastersubgrid.bf(m,k,1)<0
            pa = elemnodes(perm1(:,m),:,k);
            for i=1:size(perm1,2)
                for j=1:size(mastersubgrid.tf,1)
                    pb = boudnodes(:,:,j,i);
                    if norm(pa-pb)<1e-10
                        mastersubgrid.bf(m,k,2) = i;
                        mastersubgrid.bf(m,k,3) = j;
                        mastersubgrid.cf(j,i,1) = k;
                        mastersubgrid.cf(j,i,2) = m;
                    end
                end
            end
        end
    end
end

mastersubgrid.ibdu = zeros(npf*nsf*nfe,1);
mastersubgrid.ibduhat = mastersubgrid.ibdu;
for ls=1:nfe
    for ks=1:nsf
        k  = mastersubgrid.cf(ks,ls,1); % element
        is = mastersubgrid.cf(ks,ls,2); % sub-face
        im = (ls-1)*npf*nsf+(ks-1)*npf+1:(ls-1)*npf*nsf+ks*npf;        
        mastersubgrid.ibdu(im) = mastersubgrid.perm(:,is)+(k-1)*npv;        
        mastersubgrid.ibduhat(im) = ((is-1)*npf+1:is*npf)+(k-1)*(npf*nfe);         
        in = (ls-1)*nqf*nsf+(ks-1)*nqf+1:(ls-1)*nqf*nsf+ks*nqf;
        mastersubgrid.ibdgeom(in) = mastersubgrid.permgeom(:,is)+(k-1)*npm;
    end
end

mastersubgrid.nt  = size(mastersubgrid.t,1);      % number of elements
mastersubgrid.ngv = size(masterelem.gpvl,1);   % number of gasss points per element
mastersubgrid.ngf = size(masterelem.gpfc,1);   % number of gasss points per face
mastersubgrid.npv = size(masterelem.plocvl,1); % number of nodal points per element
mastersubgrid.npf = size(masterelem.plocfc,1); % number of nodal points per face
mastersubgrid.npm = size(mastergeom.plocvl,1); % number of geometry points per element

function dgnodes = getbounodes(p,plocal)

np=size(plocal,1);
nd=size(plocal,2);
dgnodes = zeros(np,nd+1);

% Allocate nodes
if nd==1
    xi  = plocal(:,1);
    philocal(:,1) = 1 - xi;
    philocal(:,2) = xi;   
    for i=1:np    
        dgnodes(i,:) = p(1,:)*philocal(i,1) + p(2,:)*philocal(i,2);
    end    
elseif nd==2 && npv==3 % tri
    xi  = plocal(:,1);
    eta = plocal(:,2);    
    philocal(:,1) = 1 - xi - eta;
    philocal(:,2) = xi;
    philocal(:,3) = eta;
    for i=1:np    
        dgnodes(i,:) = p(1,:)*philocal(i,1) + p(2,:)*philocal(i,2) + p(3,:)*philocal(i,3);
    end    
elseif nd==2 && npv==4 % quad
    xi  = plocal(:,1);
    eta = plocal(:,2);
    philocal(:,1) = (1-xi).*(1-eta);
    philocal(:,2) = xi.*(1-eta);
    philocal(:,3) = xi.*eta;
    philocal(:,4) = (1-xi).*eta;
    for i=1:np    
        dgnodes(i,:) = p(1,:)*philocal(i,1) + p(2,:)*philocal(i,2) + ...
                       p(3,:)*philocal(i,3) + p(4,:)*philocal(i,4);
    end    
end
    



