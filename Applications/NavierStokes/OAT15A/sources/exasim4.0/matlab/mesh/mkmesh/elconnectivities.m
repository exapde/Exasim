function [elcon,ndof,edg] = elconnectivities(mesh)
%MKELCON compute element connectivities 
%   [ELCON,NDOF] = GETELCON(MESH)
%
%      MESH:      MESH STRUCTURE
%      ELCON:     Element connectivies
%      NDOF:      Number of degrees of freedom
%--------------------------------------------------------------------------
% REVISION HISTORY:
% When     Who               What
% 28Nov12  Hemant Chaurasia  Added all(in==0) case for periodic boundaries
% 28Nov12  Hemant Chaurasia  Added porder==0 code block from HDGv3
%--------------------------------------------------------------------------

dim      = mesh.nd;

if dim == 1
    elcon = mesh.t2f';    
    return;
end

if mesh.porder==0
    npf   = 1;
    ne    = size(mesh.t,1);
    nfv   = size(mesh.t2f,2); % number of faces per element
    nf    = size(mesh.f,1);
    ncf = dim;               % number of corners of a face
    if dim == 3 && mesh.elemtype  % hex element
        ncf = dim+1;
    end    
    elcon = zeros(npf,nfv,ne);    
    for i = 1:nf        
        fe = mesh.f(i,ncf+1:ncf+2); % neighboring elements of face i        
        if1 = (mesh.t2f(fe(1),:)==i);           % location of face i on element fe(1)
        elcon(:,if1,fe(1)) = i;   % assign dof numbering of face i to elcon from element fe(1) 
        if fe(2)>0        
            if2 = (mesh.t2f(fe(2),:)==i);       % location of face i on element fe(2)
            elcon(:,if2,fe(2)) = i; % assign dof numbering of face i to elcon from element fe(2) 
        end                        
    end    
    elcon = reshape(elcon,[npf*nfv ne]);
    return;
end

porder   = mesh.porder;
elemtype = mesh.elemtype;
perm     = mesh.perm;   

ne       = size(mesh.t,1);
nfv      = size(mesh.t2f,2); % number of faces per element
[nf,nf2] = size(mesh.f); % total number of faces
ncf = dim;               % number of corners of a face
if dim == 3 && elemtype  % hex element
    ncf = dim+1;
end

[philocfc,philocvl] = localbasis(porder,dim,elemtype);
[npv,nnv] = size(philocvl);
npf       = size(philocfc,1);

% get dgnodes
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*mesh.p(mesh.t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end

elcon = zeros(npf,nfv,ne);
edg   = [];
ndof  = 0;
for i = 1:nf
    fn = mesh.f(i,1:ncf);       % corners of face i 
    fe = mesh.f(i,ncf+1:ncf+2); % neighboring elements of face i
    ff = mesh.f(i,nf2);         % 0 -> hdg face or 1 -> edg face
    
    pf = philocfc*mesh.p(fn,:); % nodal points on face i
    
    if (ff==0)                  % hdg face         
        tof = ndof + (1:npf);   % dof numbering on face i
        ndof = ndof + npf;      % number of degrees of freedom 
    elseif (ff==1)              % edg face 
        if isempty(edg)
            edg = pf;           % edg nodes on faces
            tof = ndof + (1:npf);            
            ieg = tof;          % dof numbering of edg nodes
            ndof = ndof + npf; 
        else                        
            in = xiny(pf,edg);     % find which rows of pf are in edg 
            jn = (in==0);          % new edg nodes have in==0
            kn = (in>0);           % old edg nodes have in>0  
            neg = sum(jn);         % number of new edg nodes                        
            
            edg = [edg; pf(jn,:)]; % update edg with new edg nodes            
           
            tof(jn) = ndof+(1:neg); % tof for new edge nodes            
            tof(kn) = ieg(in(kn));  % tof for old edge nodes                                                      
            
            ieg = [ieg ndof+(1:neg)]; % update ieg with dof numbering of new edge nodes 
            ndof = ndof + neg;        % update ndof     
        end
    end
        
    if1 = (mesh.t2f(fe(1),:)==i);           % location of face i on element fe(1)
    dg1 = squeeze(dgnodes(perm(:,if1),:,fe(1)));  % dg nodes on face i from element fe(1)
    in = xiny(dg1,pf);              % match pf to dg1  
    try
        elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1) 
    catch err
        1;
    end
        
    if fe(2)>0
        if2 = (mesh.t2f(fe(2),:)==i);       % location of face i on element fe(2)
        dg2 = squeeze(dgnodes(perm(:,if2),:,fe(2))); % dg nodes on face i from element fe(2)
        in = xiny(dg2,pf);            % match pf to dg2  
        if all(in==0)   % case for faces on periodic boundaries (Hemant Chaurasia)
            in1 = xiny(dg1,pf);
            in = in1(size(dg1,1):-1:1);  % reverse of in1
        end
        elcon(:,if2,fe(2)) = tof(in); % assign dof numbering of face i to elcon from element fe(2) 
    end                
end

elcon = reshape(elcon,[npf*nfv ne]);


function [philocfc,philocvl] = localbasis(porder,dim,elemtype) 

[plocvl,tlocal,plocfc] = masternodes(porder,dim,elemtype,0);

if dim==2 && elemtype==0      % tri
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = 1 - xi - eta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
elseif dim==2 && elemtype==1  % quad
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = (1-xi).*(1-eta);
    philocvl(:,2) = xi.*(1-eta);
    philocvl(:,3) = xi.*eta;
    philocvl(:,4) = (1-xi).*eta;
elseif dim==3 && elemtype==0  % tet
    xi  = plocfc(:,1);
    eta = plocfc(:,2);    
    philocfc(:,1) = 1 - xi - eta;
    philocfc(:,2) = xi;
    philocfc(:,3) = eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = 1 - xi - eta - zeta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
    philocvl(:,4) = zeta;
elseif dim==3 && elemtype==1   % hex
    xi  = plocfc(:,1);
    eta = plocfc(:,2);
    philocfc(:,1) = (1-xi).*(1-eta);
    philocfc(:,2) = xi.*(1-eta);
    philocfc(:,3) = xi.*eta;
    philocfc(:,4) = (1-xi).*eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocvl(:,2) = xi.*(1-eta).*(1-zeta);
    philocvl(:,3) = xi.*eta.*(1-zeta);
    philocvl(:,4) = (1-xi).*eta.*(1-zeta);    
    philocvl(:,5) = (1-xi).*(1-eta).*(zeta);
    philocvl(:,6) = xi.*(1-eta).*(zeta);
    philocvl(:,7) = xi.*eta.*(zeta);
    philocvl(:,8) = (1-xi).*eta.*(zeta);        
end


function in = xiny(x,y)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

[m,dim] = size(x);

in = zeros(m,1);
% if dim==2
%     for j=1:m
%         d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
%         [md,id] = min(d2);
%         if md<1e-12, in(j)=id; end
%     end
% else
%     for j=1:m
%         d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
%         [md,id] = min(d2);
%         if md<1e-12, in(j)=id; end
%     end
% end

if dim == 2
    y1Rep = repmat(y(:,1),[1,m]);
    y2Rep = repmat(y(:,2),[1,m]);
    x1Rep = repmat(x(:,1)',[size(y,1),1]);
    x2Rep = repmat(x(:,2)',[size(y,1),1]);

    dis = (y1Rep-x1Rep).^2 + (y2Rep-x2Rep).^2;
    [md,id] = min(dis,[],1);
    aux = md < (1e-10)^2;
    in(aux) = id(aux);
else
    y1Rep = repmat(y(:,1),[1,m]);
    y2Rep = repmat(y(:,2),[1,m]);
    y3Rep = repmat(y(:,3),[1,m]);
    x1Rep = repmat(x(:,1)',[size(y,1),1]);
    x2Rep = repmat(x(:,2)',[size(y,1),1]);
    x3Rep = repmat(x(:,3)',[size(y,1),1]);

    dis = (y1Rep-x1Rep).^2 + (y2Rep-x2Rep).^2 + (y3Rep-x3Rep).^2;
    [md,id] = min(dis,[],1);
    aux = md < (1e-10)^2;
    in(aux) = id(aux);
end
