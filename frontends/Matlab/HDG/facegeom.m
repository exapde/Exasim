function [pg, nlg, jac] = facegeom(shapgeomft,dgnodes,perm)
% FACEGEOM computes dg nodes, Jacobian determinant and normal vectors at Gauss points  

%   [pg, nlg, jac] = facegeom(shapgeomft,dgnodes,perm)
%
%    SHAPGEOMFT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PERM       :  Indices of the boundary nodes
%    PG         :  Physical nodes at Gauss points 
%    nlg        :  Normal vector at Gauss points
%    jac        :  Determinant of the Jacobian mapping 

ne    = size(dgnodes,2);
nq    = size(dgnodes,3);
ngf   = size(shapgeomft,1);
npf   = size(shapgeomft,2);
nd    = size(shapgeomft,3);
nfe   = size(perm,2);

perm  = perm(:,:,1);
pn = reshape(dgnodes(perm,:,:),[npf nfe*ne*nq]);

if nd>1
    dshapft  = reshape(permute(shapgeomft(:,:,2:nd),[1 3 2]),[ngf*(nd-1) npf]);
    dpg = dshapft*pn(:,1:nfe*ne*nd);
    dpg = permute(reshape(dpg,[ngf nd-1 nfe ne nd]), [1 3 4 5 2]);    
    dpg = reshape(dpg,[ngf*nfe*ne,nd,nd-1]);    
end

shapgeomft   = shapgeomft(:,:,1);
pg = shapgeomft*pn;
pg = reshape(pg,[ngf*nfe*ne nq]);

switch nd
    case 1
        jac = [ones(ne,1); ones(ne,1)];
        jac = jac(:);
        nlg = [-ones(1,ne); ones(1,ne)];        
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


