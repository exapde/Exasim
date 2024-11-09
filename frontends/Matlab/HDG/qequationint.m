function [MinvC, MinvE, Mass, Minv] = qequationint(master, dgnodes)

% get dimensions
nd    = master.nd;
npe   = master.npe;
nge   = master.nge;
ngf   = master.ngf;
npf   = master.npf;
nfe   = size(master.perm,2);
ne    = size(dgnodes,3);

% get shape functions and their derivatives
perm = master.perm;
dshapvt = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[nge*nd npe]);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npe*npe nge (nd+1)]);
dshapft = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

Mass = zeros(npe, npe, ne);
Minv = Mass;
MinvC = zeros(npe, npe, ne, nd);
MinvE = zeros(npe, npf*nfe, ne, nd);

% loop through each element        
for k=1:1:ne 
    % get dg nodes
    dg = dgnodes(:,:,k);    
    
    % get dg nodes on faces
    pn = reshape(dg(perm,1:nd,:),[npf nfe*nd]);
    
    % compute volumetic Jacobian matrix 
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[nge nd nd]);        
    % compute the determintant and inverse
    [jac,Xx] = volgeom(Jg);    
            
    dpg = dshapft*pn;
    dpg = permute(reshape(dpg,[ngf nd-1 nfe nd]), [1 3 4 2]);    
    dpg = reshape(dpg,[ngf*nfe,nd,nd-1]); 
    % get the normal vectors and Jacobian determinants on faces
    [nlg, jacf] = facegeom(dpg,nd);
    
    % mass matrix
    M = reshape(shapvgdotshapvl(:,:,1)*jac,[npe npe]);    
    Mi= inv(M);
    
    Mass(:,:,k) = M;
    Minv(:,:,k) = Mi;
    
    % convection matrices
    C = zeros(npe,npe,nd); 
    for i=1:nd
        C(:,:,i) = reshape(shapvgdotshapvl(:,:,2)*Xx(:,i,1),[npe npe]);        
        for j=2:nd
            C(:,:,i) = C(:,:,i) + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,i,j),[npe npe]);
        end
        MinvC(:,:,k,i) = -Mi*C(:,:,i);
    end                
        
    % face matrices
    E = zeros(npe,npf*nfe,nd);    
    for i=1:nd    
        njc = reshape(nlg(:,i).*(jacf),[ngf,nfe]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe]);   
        for j=1:nfe
            E(perm(:,j),(j-1)*npf+1:j*npf,i) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end            
        MinvE(:,:,k,i) = Mi*E(:,:,i);
    end        
end

function [jac,Xx] = volgeom(Jg)

ngv = size(Jg,1); 
nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;
        Xx = -ones(ngv,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = -Jg(:,2,2);
        Xx(:,2,1) = Jg(:,2,1);
        Xx(:,1,2) = Jg(:,1,2);
        Xx(:,2,2) = -Jg(:,1,1);
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,3).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,3);
        Xx(:,2,1) = Jg(:,2,1).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,1);
        Xx(:,3,1) = Jg(:,2,2).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,2);
        Xx(:,1,2) = Jg(:,1,2).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,2);
        Xx(:,2,2) = Jg(:,1,3).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,3);
        Xx(:,3,2) = Jg(:,1,1).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,1);
        Xx(:,1,3) = Jg(:,1,3).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,3);
        Xx(:,2,3) = Jg(:,1,1).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,1);
        Xx(:,3,3) = Jg(:,1,2).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,2);
    otherwise
        error('Dimension is not implemented');
end


function [nlg, jac] = facegeom(dpg,nd)

switch nd
    case 1
        jac = 1;
        nlg = [1; -1];        
        nlg = nlg(:);
    case 2
        jac = sqrt(dpg(:,1).^2+dpg(:,2).^2);
        nlg   = [dpg(:,2),-dpg(:,1)];
        nlg   = bsxfun(@rdivide, nlg, jac);
    case 3
        nlg(:,1) = dpg(:,2,1).*dpg(:,3,2) - dpg(:,3,1).*dpg(:,2,2);
        nlg(:,2) = dpg(:,3,1).*dpg(:,1,2) - dpg(:,1,1).*dpg(:,3,2);
        nlg(:,3) = dpg(:,1,1).*dpg(:,2,2) - dpg(:,2,1).*dpg(:,1,2);
        jac = sqrt(nlg(:,1).^2+nlg(:,2).^2+nlg(:,3).^2);
        nlg   = bsxfun(@rdivide, nlg, jac);
    otherwise
        error('Dimension is not implemented');
end








