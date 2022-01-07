function master=mkmasterelement(dim,porder,morder,pgauss,pgaussR,elemtype,nodetype)
%MKMASTER  Create master element structure
%    MASTER=MKMASTER(MESH)
%
%      MESH:      Mesh data structure
%      PGAUSS:    Degree of the polynomila to be integrated exactly
%                 (default: PGAUSS = 2*MESH.PORDER)
%

if isempty(morder)
    morder = porder;
end
if isempty(pgauss)
    pgauss = max(2*porder,1); 
end
if isempty(pgaussR)
    pgaussR = max(2*porder,1); 
end

master.nd     = dim;      % problem dimension
master.porder = porder;   % polynomial degree
master.morder = morder;   % polynomial degree
master.pgauss = pgauss;
master.pgaussR = pgaussR;
master.elemtype = elemtype;
master.nodetype = nodetype;

% vertex basic functions at the DG nodes
[master.philocvl,master.philocfc] = localbasis(porder,dim,elemtype); 

% node positions on the master element
[master.plocvl,master.tlocvl,master.plocfc,master.tlocfc,permnode,permedge,permface,master.face] = mkmasternodes(porder,dim,elemtype,nodetype);
master.permnode = permnode; % positions of nodes at the vertices
master.permedge = permedge; % positions of nodes on the edges
master.permface = permface; % positions of nodes on the faces
% positions of nodes on the faces
if dim==1
    master.perm = permnode;
elseif dim==2
    master.perm = permedge;
elseif dim==3
    master.perm = permface;
end
if iscell(master.perm); master.perm = cell2mat(master.perm); end

% geometry nodes on the master element and face
[master.plocmv,~,master.plocmf,~,~,~,master.permgeom] = mkmasternodes(morder,dim,elemtype,nodetype);

% Gauss points and weights on the master volume element  
[master.gpvl,master.gwvl] = gaussquad(pgauss,dim,elemtype);
[master.gpvlR,master.gwvlR] = gaussquad(pgaussR,dim,elemtype);

% master.gpnvl = master.plocvl;
% master.gwnvl = (1/size(master.plocvl,1)) * ones(size(master.plocvl,1),1);

% shape functions and derivatives on the master volume element  
master.shapvl = mkshape(porder,master.plocvl,master.gpvl,elemtype);
master.shapnl = mkshape(porder,master.plocvl,master.plocvl,elemtype);
master.shapvlR = mkshape(porder,master.plocvl,master.gpvlR,elemtype);
% master.shapnvl = mkshape(porder,master.plocvl,master.plocvl,elemtype);

master.ploc1d = mkmasternodes(morder,1,elemtype,nodetype);
[master.gp1d,master.gw1d] = gaussquad(pgauss,1,elemtype);
master.shap1dg = mkshape(porder,master.ploc1d,master.gp1d,elemtype);
master.shap1dn = mkshape(porder,master.ploc1d,master.ploc1d,elemtype);
for d=1:2
    master.shap1dgt(:,:,d) = master.shap1dg(:,:,d)';    
    master.shap1dgw(:,:,d) = master.shap1dg(:,:,d)*diag(master.gw1d);    
    master.shap1dnt(:,:,d) = master.shap1dn(:,:,d)';   
    master.shap1dnl(:,:,d) = master.shap1dn(:,:,d);   
end

% geometry shape functions
master.shapmv = mkshape(morder,master.plocmv,master.gpvl,elemtype);
master.shapmvR = mkshape(morder,master.plocmv,master.gpvlR,elemtype);
master.shapmv = permute(master.shapmv,[2 1 3]);
master.shapmvR = permute(master.shapmvR,[2 1 3]);

ngv = size(master.gpvl,1);   % number of gasss points per element
npv = size(master.plocvl,1); % number of nodal points per element
% volume shape functions and their derivatives 
master.shapvgdotshapvl  = zeros(npv*npv,ngv,dim+1);      
for d=1:dim+1
    master.shapvt(:,:,d) = master.shapvl(:,:,d)';   
    master.shapnt(:,:,d) = master.shapnl(:,:,d)';   
    master.shapvg(:,:,d) = master.shapvl(:,:,d)*diag(master.gwvl);    
    master.shapvtR(:,:,d) = master.shapvlR(:,:,d)';
    master.shapvgR(:,:,d) = master.shapvlR(:,:,d)*diag(master.gwvlR);    
    for ii=1:npv
        for jj = 1:npv
            master.shapvgdotshapvl((ii-1)*npv+jj,:,d) = master.shapvg(jj,:,d).*master.shapvl(ii,:,1);                    
        end
    end            
end

% master.shapvndotshapvl  = zeros(npv*npv,npv,dim+1);      
% for d=1:dim+1
%     master.shapvnt(:,:,d) = master.shapnvl(:,:,d)';
%     master.shapvn(:,:,d) = master.shapnvl(:,:,d)*diag(master.gwnvl);    
%     for ii=1:npv
%         for jj = 1:npv
%             master.shapvndotshapvl((ii-1)*npv+jj,:,d) = master.shapvn(jj,:,d).*master.shapnvl(ii,:,1);                    
%         end
%     end
% end

% shape functions and derivatives on the master face element  
nfe = length(master.plocfc);
[master.gpfc,master.gwfc,master.shapfc] = gaussshapeface(master.plocfc,porder,pgauss,dim,elemtype);
[master.gpfcR,master.gwfcR,master.shapfcR] = gaussshapeface(master.plocfc,porder,pgaussR,dim,elemtype);

% geometry shape functions and derivatives on the master face element  
[~,~,master.shapmf] = gaussshapeface(master.plocmf,morder,pgauss,dim,elemtype);
[~,~,master.shapmfR] = gaussshapeface(master.plocmf,morder,pgaussR,dim,elemtype);

% face shape functions and their derivatives 
for i = 1:nfe
    master.shapfn{i} = mkshape(porder(1),master.plocfc{i},master.plocfc{i},elemtype);
    npf = size(master.shapfc{i},1);    
    ngf = size(master.shapfc{i},2);
    master.shapmf{i} = permute(master.shapmf{i},[2 1 3]);        
    master.shapmfR{i} = permute(master.shapmfR{i},[2 1 3]);        
    master.shapfgdotshapfc{i}  = zeros(npf*npf,ngf,dim);      
    for d=1:dim
        master.shapfnt{i}(:,:,d) = master.shapfn{i}(:,:,d)';
        master.shapft{i}(:,:,d) = master.shapfc{i}(:,:,d)';
        master.shapfg{i}(:,:,d) = master.shapfc{i}(:,:,d)*diag(master.gwfc{i});
        master.shapftR{i}(:,:,d) = master.shapfcR{i}(:,:,d)';
        master.shapfgR{i}(:,:,d) = master.shapfcR{i}(:,:,d)*diag(master.gwfcR{i});        
        for ii=1:npf
            for jj = 1:npf
                master.shapfgdotshapfc{i}((ii-1)*npf+jj,:,d) = master.shapfg{i}(jj,:,d).*master.shapfc{i}(ii,:,1);                    
            end
        end            
    end
end

% mass and convection matrices on the master element
master.mass = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,1))';
for ii=1:dim
    master.conv(:,:,ii) = squeeze(master.shapvl(:,:,1))*diag(master.gwvl)*squeeze(master.shapvl(:,:,ii+1))';
end   

master.npv = size(master.plocvl,1); % number of nodal points per element
master.nmv = size(master.plocmv,1); % number of nodal points per element
master.ngv = size(master.gpvl,1);   % number of gasss points per element
master.ngvR = size(master.gpvlR,1);   % number of gasss points per element
for i = 1:nfe
    master.npf(i) = size(master.plocfc{i},1); % number of nodal points per face
    master.nmf(i) = size(master.plocmf{i},1); % number of nodal points per face
    master.ngf(i) = size(master.gpfc{i},1);   % number of gasss points per face
    master.ngfR(i) = size(master.gpfcR{i},1);   % number of gasss points per face     
    master.nvf(i) = length(master.face{i}); % number of vertices per face
end

nle = size(master.permedge,2);
master.npl = zeros(nle,1);
for i = 1:nle
    master.npl(i) = length(find(master.permedge(:,i)>0));
end

function [gpfc,gwfc,shapfc] = gaussshapeface(plocfc,porder,pgauss,dim,elemtype)

nfe = length(plocfc);
switch dim
    case 1
        for i=1:nfe            
            master.gpfc{i}   = 0;
            master.gwfc{i}   = 1;
            master.shapfc{i} = 1;
        end        
    case 2 % 2D
        if elemtype==0     % tri        
            for i = 1:nfe
                % Gauss points and weights on the master face element  
                [gpfc{i},gwfc{i}] = gaussquad(pgauss(1),dim-1,elemtype);

                % shape functions and derivatives on the master face element  
                shapfc{i} = mkshape(porder(1),plocfc{i},gpfc{i},elemtype);
            end                
        elseif elemtype==1 % quad
            if length(porder)==1
                porder = porder*ones(dim,1);
            end            
            if length(pgauss)==1
                pgauss = pgauss*ones(dim,1);
            end            
            [gpfc{1},gwfc{1}] = gaussquad(pgauss(1),dim-1,elemtype);
            shapfc{1} = mkshape(porder(1),plocfc{1},gpfc{1},elemtype);

            [gpfc{2},gwfc{2}] = gaussquad(pgauss(2),dim-1,elemtype);
            shapfc{2} = mkshape(porder(2),plocfc{2},gpfc{2},elemtype);

            [gpfc{3},gwfc{3}] = gaussquad(pgauss(1),dim-1,elemtype);
            shapfc{3} = mkshape(porder(1),plocfc{3},gpfc{3},elemtype);

            [gpfc{4},gwfc{4}] = gaussquad(pgauss(2),dim-1,elemtype);
            shapfc{4} = mkshape(porder(2),plocfc{4},gpfc{4},elemtype);
        end        
    case 3 %3D
        if elemtype==0     % tet
            for i = 1:nfe
                % Gauss points and weights on the master face element  
                [gpfc{i},gwfc{i}] = gaussquad(pgauss(1),dim-1,elemtype);

                % shape functions and derivatives on the master face element  
                shapfc{i} = mkshape(porder(1),plocfc{i},gpfc{i},elemtype);
            end                
        elseif elemtype==1 % hex
            if length(porder)==1
                porder = porder*ones(dim,1);
            end            
            if length(pgauss)==1
                pgauss = pgauss*ones(dim,1);
            end                        
            [gpfc{1},gwfc{1}] = gaussquad([pgauss(2) pgauss(1)],dim-1,elemtype);
            shapfc{1} = mkshape([porder(2) porder(1)],plocfc{1},gpfc{1},elemtype);

            [gpfc{2},gwfc{2}] = gaussquad([pgauss(1) pgauss(2)],dim-1,elemtype);
            shapfc{2} = mkshape([porder(1) porder(2)],plocfc{2},gpfc{2},elemtype);

            [gpfc{3},gwfc{3}] = gaussquad([pgauss(1) pgauss(3)],dim-1,elemtype);
            shapfc{3} = mkshape([porder(1) porder(3)],plocfc{3},gpfc{3},elemtype);

            [gpfc{4},gwfc{4}] = gaussquad([pgauss(1) pgauss(3)],dim-1,elemtype);
            shapfc{4} = mkshape([porder(1) porder(3)],plocfc{4},gpfc{4},elemtype);

            [gpfc{5},gwfc{5}] = gaussquad([pgauss(2) pgauss(3)],dim-1,elemtype);
            shapfc{5} = mkshape([porder(2) porder(3)],plocfc{5},gpfc{5},elemtype);

            [gpfc{6},gwfc{6}] = gaussquad([pgauss(2) pgauss(3)],dim-1,elemtype);
            shapfc{6} = mkshape([porder(2) porder(3)],plocfc{6},gpfc{6},elemtype);
        elseif elemtype==2 % prism
            if length(porder)==1
                porder = [porder porder];
            end                        
            if length(pgauss)==1
                pgauss = [pgauss pgauss];
            end                                                
            [gpfc{1},gwfc{1}] = gaussquad(pgauss(1),dim-1,0);
            shapfc{1} = mkshape(porder(1),plocfc{1},gpfc{1},0);

            [gpfc{2},gwfc{2}] = gaussquad(pgauss(1),dim-1,0);
            shapfc{2} = mkshape(porder(1),plocfc{2},gpfc{2},0);

            [gpfc{3},gwfc{3}] = gaussquad([pgauss(1) pgauss(2)],dim-1,1);
            shapfc{3} = mkshape([porder(1) porder(2)],plocfc{3},gpfc{3},1);

            [gpfc{4},gwfc{4}] = gaussquad([pgauss(1) pgauss(2)],dim-1,1);
            shapfc{4} = mkshape([porder(1) porder(2)],plocfc{4},gpfc{4},1);

            [gpfc{5},gwfc{5}] = gaussquad([pgauss(1) pgauss(2)],dim-1,1);
            shapfc{5} = mkshape([porder(1) porder(2)],plocfc{5},gpfc{5},1);
        elseif elemtype==3 % pyramid
            [gpfc{1},gwfc{1}] = gaussquad(pgauss,dim-1,1);
            shapfc{1} = mkshape(porder,plocfc{1},gpfc{1},1);
            for i = 2:nfe
                [gpfc{i},gwfc{i}] = gaussquad(pgauss,dim-1,0);
                shapfc{i} = mkshape(porder,plocfc{i},gpfc{i},0);
            end
        end
    otherwise
        error('Only can handle dim=1, dim=2 or dim=3');
end            
