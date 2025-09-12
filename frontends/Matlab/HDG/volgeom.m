function [pg, Xx, jac] = volgeom(shapgeomvt,dgnodes)
% VOLGEOM computes physical nodes, Jacobian determinant and inverse at Gauss points  

%   [pg, Xx, jac] = volgeom(shapgeomvt,dgnodes)
%
%    SHAPGEOMVT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PG         :  Physical nodes at Gauss points 
%    Xx         :  Minus the inverse of the Jacobian mapping times the determinant
%    jac        :  Determinant of the Jacobian mapping 

nge   = size(shapgeomvt,1);   % number of gauss points
npe   = size(shapgeomvt,2);
nd    = size(shapgeomvt,3)-1;
nq    = size(dgnodes,3);
ne    = size(dgnodes,2);

dshapvt   = reshape(permute(shapgeomvt(:,:,2:nd+1),[1 3 2]),[nge*nd npe]);
shapgeomvt    = shapgeomvt(:,:,1); 

% compute dg nodes at Gauss points: x = sum_i phi_i(xi) *x_i
pn = reshape(dgnodes,[npe ne*nq]);
pg = shapgeomvt*pn;
pg = reshape(pg,[nge*ne nq]);

% compute the Jacobian matrix at Gauss points: dx/dxi
Jg = dshapvt*pn(:,1:ne*nd);
Jg = permute(reshape(Jg,[nge nd ne nd]),[1 3 2 4]);
Jg = reshape(Jg,[nge,ne,nd,nd]);

% jac: the determinant of the Jacobian matrix at Gauss points
% Xx = inv(Jg)*jac: the inverse of the Jacobian matrix times 
%                   the determinant of the Jacobian matrix 
switch nd
    case 1
        jac = Jg;
        Xx = ones(nge,ne);
    case 2
        jac = Jg(:,:,1,1).*Jg(:,:,2,2) - Jg(:,:,1,2).*Jg(:,:,2,1);
        Xx(:,:,1,1) = Jg(:,:,2,2);
        Xx(:,:,2,1) = -Jg(:,:,2,1);
        Xx(:,:,1,2) = -Jg(:,:,1,2);
        Xx(:,:,2,2) = Jg(:,:,1,1);
    case 3
        jac = Jg(:,:,1,1).*Jg(:,:,2,2).*Jg(:,:,3,3) - Jg(:,:,1,1).*Jg(:,:,3,2).*Jg(:,:,2,3)+ ...
              Jg(:,:,2,1).*Jg(:,:,3,2).*Jg(:,:,1,3) - Jg(:,:,2,1).*Jg(:,:,1,2).*Jg(:,:,3,3)+ ...
              Jg(:,:,3,1).*Jg(:,:,1,2).*Jg(:,:,2,3) - Jg(:,:,3,1).*Jg(:,:,2,2).*Jg(:,:,1,3);            
        Xx(:,:,1,1) = Jg(:,:,2,2).*Jg(:,:,3,3) - Jg(:,:,2,3).*Jg(:,:,3,2);
        Xx(:,:,2,1) = Jg(:,:,2,3).*Jg(:,:,3,1) - Jg(:,:,2,1).*Jg(:,:,3,3);
        Xx(:,:,3,1) = Jg(:,:,2,1).*Jg(:,:,3,2) - Jg(:,:,2,2).*Jg(:,:,3,1);
        Xx(:,:,1,2) = Jg(:,:,1,3).*Jg(:,:,3,2) - Jg(:,:,1,2).*Jg(:,:,3,3);
        Xx(:,:,2,2) = Jg(:,:,1,1).*Jg(:,:,3,3) - Jg(:,:,1,3).*Jg(:,:,3,1);
        Xx(:,:,3,2) = Jg(:,:,1,2).*Jg(:,:,3,1) - Jg(:,:,1,1).*Jg(:,:,3,2);
        Xx(:,:,1,3) = Jg(:,:,1,2).*Jg(:,:,2,3) - Jg(:,:,1,3).*Jg(:,:,2,2);
        Xx(:,:,2,3) = Jg(:,:,1,3).*Jg(:,:,2,1) - Jg(:,:,1,1).*Jg(:,:,2,3);
        Xx(:,:,3,3) = Jg(:,:,1,1).*Jg(:,:,2,2) - Jg(:,:,1,2).*Jg(:,:,2,1);
    otherwise
        error('Dimension is not implemented');
end

jac = reshape(jac, [nge*ne 1]);
