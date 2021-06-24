function QDG = getq3(master, mesh, UDG, SH, fc_q)

if nargin < 5; SH = []; end
if nargin < 6; fc_q = 1; end

% get dimensions
nd    = master.nd;
npv   = master.npv;
ngv   = master.ngv;
ngf   = master.ngf;
npf   = master.npf;
%npe   = master.npe;
nfe   = size(master.perm,2);
ne    = size(UDG,3);
nch    = size(UDG,2);
%nch   = size(UH,2);
ncq   = nch*nd;

%UH    = UH(:,mesh.elcon);
% UH = reshape(permute(UH,[2 1 3]),[nch npf*mesh.nf]);
% UH = permute(reshape(UH(:,mesh.elcon),[nch nfe*npf ne]),[2 1 3]);

master.shapfgdotshapfc = zeros(npf*npf,ngf,nd);
for d=1:nd
    master.shapft(:,:,d) = master.shapfc(:,:,d)';
    master.shapfg(:,:,d) = master.shapfc(:,:,d)*diag(master.gwfc);
    for ii=1:npf
        for jj = 1:npf
            master.shapfgdotshapfc((ii-1)*npf+jj,:,d) = master.shapfg(jj,:,d).*master.shapfc(ii,:,1);
        end
    end
end

% get shape functions and their derivatives
perm = master.perm;
dshapvt = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
shapvgdotshapvl = reshape(master.shapvgdotshapvl,[npv*npv ngv (nd+1)]);
dshapft = reshape(permute(master.shapft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
shapfgdotshapfc = master.shapfgdotshapfc(:,:,1);

% allocate memory
QDG = zeros(npv,ncq,ne);

% loop through each element        
for k=1:1:ne 
    % get dg nodes
    dg = mesh.dgnodes(:,:,k);    
    
    % get dg nodes on faces
    pn = reshape(dg(perm,:,:),[npf nfe*nd]);
    
    % compute volumetic Jacobian matrix 
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    % compute the determintant and inverse
    [jac,Xx] = volgeom(Jg);
        
    dpg = dshapft*pn;
    dpg = permute(reshape(dpg,[ngf nd-1 nfe nd]), [1 3 4 2]);    
    dpg = reshape(dpg,[ngf*nfe,nd,nd-1]); 
    % get the normal vectors and Jacobian determinants on faces
    [nlg, jacf] = facegeom(dpg,nd);
    
    % mass matrix
    M = reshape(shapvgdotshapvl(:,:,1)*jac,[npv npv]);    
    Mi= inv(M)/fc_q;
    
    % convection matrices
    C = zeros(npv,npv,nd); 
    for i=1:nd
        C(:,:,i) = reshape(shapvgdotshapvl(:,:,2)*Xx(:,i,1),[npv npv]);        
        for j=2:nd
            C(:,:,i) = C(:,:,i) + reshape(shapvgdotshapvl(:,:,j+1)*Xx(:,i,j),[npv npv]);
        end            
    end    
    C = reshape(permute(C,[1 3 2]),[npv*nd npv]);
    
    % face matrices
    E = zeros(npv,npf*nfe,nd);    
    for i=1:nd    
        njc = reshape(nlg(:,i).*(-jacf),[ngf,nfe]);
        wrk = reshape(shapfgdotshapfc*njc,[npf npf*nfe]);   
        for j=1:nfe
            E(perm(:,j),(j-1)*npf+1:j*npf,i) = wrk(1:npf,(j-1)*npf+1:j*npf,:);      
        end            
    end    
    E = reshape(permute(E,[1 3 2]),[npv*nd npf*nfe]);
        
    % obtain the current solution and numerical trace    
    u    = UDG(:,1:nch,k);    
    %uh   = UH(:,:,k);         
    %sh   = SH(:,end-ncq+1:end,k);
    uh = zeros(npf,nch,nfe);    
    for i=1:nfe % for each face i of element k       
        fki = mesh.t2f(k,i); % global face fki       
        eki = mesh.f(fki,:); % 
        e1 = eki(end-1); % element e1 of face fki
        e2 = eki(end);   % element e2 of face fki
        if e2>0 % if it is an interior face
            if (e1==k) % element k coincides with e1
                u1 = u(perm(:,i),:);
                x1 = mesh.dgnodes(perm(:,i),1:nd,k); 
                j = (mesh.t2f(e2,:)==fki); % location of face fki on element e2                                
                x2 = mesh.dgnodes(perm(:,j),1:nd,e2); 
                in = xiny(round(x1,6),round(x2,6));
                u2 = UDG(perm(in,j),1:nch,e2);                
            elseif (e2==k) % element k coincides with e2
                u1 = u(perm(:,i),:);
                x1 = mesh.dgnodes(perm(:,i),1:nd,k);
                j = (mesh.t2f(e1,:)==fki); % location of face fki on element e1
                x2 = mesh.dgnodes(perm(:,j),1:nd,e1); 
                in = xiny(round(x1,6),round(x2,6));
                u2 = UDG(perm(in,j),1:nch,e1);
            else
                error('something wrong');
            end
            uh(:,:,i) = (u1+u2)/2;        
        else
            u1 = u(perm(:,i),:);            
            %x1 = mesh.dgnodes(perm(:,i),1:nd,k); 
            %uh(:,:,i) = boundaryconditions(u1, x1);
            uh(:,:,i) = u1;
        end        
    end
    uh = reshape(permute(uh,[1 3 2]),[npf*nfe nch]);
    
    % compute dq        
    Cu   = reshape(C*u,[npv nd*nch]);
    Euh  = reshape(E*uh,[npv nd*nch]);
    MiCu = permute(reshape(Mi*Cu,[npv nd nch]), [1 3 2]);
    MiEu = permute(reshape(Mi*Euh,[npv nd nch]), [1 3 2]);
    q    = reshape(MiEu,[npv nch*nd]) - reshape(MiCu,[npv nch*nd]);         
    QDG(:,:,k) = q; 
end

if isempty(SH)==0 && fc_q ~= 0
    QDG = QDG + SH(:,end-ncq+1:end,:)/fc_q;
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







