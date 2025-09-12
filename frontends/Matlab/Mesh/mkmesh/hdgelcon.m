function [elcon,ndof,hdg] = hdgelcon(mesh)

dim      = mesh.nd;
porder   = mesh.porder;
elemtype = mesh.elemtype;
perm = mesh.perm(:);

ne       = size(mesh.t,1);   % number of elements
nfe      = size(mesh.t2f,2); % number of faces per element
nf       = size(mesh.f,1);   % number of faces
ncf = dim;                   % number of corners of a face
if dim == 3 && elemtype      % hex element
    ncf = dim+1;
end

% get local polymomials on element and face
[philocfc,philocvl] = localbasis(porder,dim,elemtype);
[npv,nnv] = size(philocvl);   % number of nodes on an element and number of vertices of an element
npf       = size(philocfc,1); % number of nodes on a face

% get dg nodes
dgnodes = zeros(npv,dim,ne);
for d=1:dim
  for n=1:nnv
    dp=philocvl(:,n)*mesh.p(mesh.t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end
snap = 1e-8;
dgnodes = round(dgnodes/snap)*snap;
%dgnodes = reshape(permute(dgnodes(perm,1:dim,:),[1,3,2]),[npf*nfe*ne,dim]);
dgnodes = reshape(dgnodes(perm,1:dim,:),[npf,nfe,dim,ne]);

% get hdg nodes
pf = reshape(mesh.p(mesh.f(:,1:ncf)',:),[ncf nf*dim]); % nodal points on faces    
hdg = permute(reshape(philocfc*pf,[npf nf dim]),[1 3 2]);
hdg = round(hdg/snap)*snap;
ndof = npf*nfe;

% get element-to-entity connectivities 
elcon = zeros(npf,nfe,ne);
for i = 1:nf
    fe = mesh.f(i,ncf+1:ncf+2); % neighboring elements of face i        
    pf = hdg(:,:,i);
    
    if1 = mesh.t2f(fe(1),:)==i;           % location of face i on element fe(1)
    dg1 = reshape(dgnodes(:,if1,:,fe(1)),[npf dim]); 
    
    tof = ((i-1)*npf+1):1:(i*npf);
    in = xiny(dg1,pf);              % match pf to dg1          
    elcon(:,if1,fe(1)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)             
    
    if fe(2)>0
        if2 = find(mesh.t2f(fe(2),:)==i);           % location of face i on element fe(2)
        dg2 = reshape(dgnodes(:,if2,:,fe(2)),[npf dim]); 
        in = xiny(dg2,pf);              % match pf to dg1          
        elcon(:,if2,fe(2)) = tof(in);   % assign dof numbering of face i to elcon from element fe(1)                     
    end
end
elcon = reshape(elcon,[npf*nfe ne]);

function [philocfc,philocvl] = localbasis(porder,dim,elemtype) 

[plocvl,~,plocfc] = mkmasternodes(porder,dim,elemtype,0);

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

