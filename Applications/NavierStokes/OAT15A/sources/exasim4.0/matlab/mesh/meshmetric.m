function [M,Minv] = meshmetric(dgnodes,porder,elemtype,nodetype)

[np,nd,ne] = size(dgnodes);
nm = np*ne;

plocvl = mkmasternodes(porder,nd,elemtype,nodetype);
shapnl = mkshape(porder,plocvl,plocvl,elemtype);
shapnt = permute(shapnl,[2 1 3]);

[J,tm] = voljac(shapnt,dgnodes);
if max(abs(tm(:)-dgnodes(:)))>1e-12
    error('something wrong');
end

A = eye(nd,nd);
if (nd==2) && (elemtype==0)
    A(1,2) = -1.0/sqrt(3.0);
    A(2,2) = 2.0/sqrt(3.0);
elseif (nd==3) && (elemtype==0)
    A(1,2) = -1.0/sqrt(3.0);
    A(2,2) = sqrt(4.0/3.0);
    A(1,3) = -sqrt(1.0/6.0);
    A(2,3) = -sqrt(1.0/6.0);
    A(3,3) = sqrt(3.0/2.0);
end

M = zeros(nd,nd,nm);
Minv = M;
for i = 1:nm
    B = J(:,:,i)*A;
    M(:,:,i) = B'*B;
    Minv(:,:,i) = inv(M(:,:,i));
end
M = reshape(M,[nd*nd np ne]); 
M = permute(M,[2 1 3]);
Minv = reshape(Minv,[nd*nd np ne]); 
Minv = permute(Minv,[2 1 3]);


function [Jg, pg] = voljac(shapgeomvt,dgnodes)
% VOLGEOM computes physical nodes, Jacobian determinant and inverse at Gauss points  

%   [pg, Xx, jac] = volgeom(shapgeomvt,dgnodes)
%
%    SHAPGEOMVT :  Shape functions and derivatives at Gauss points
%    DGNODES    :  Geometry DG nodes 
%    PG         :  Physical nodes at Gauss points 
%    Xx         :  Minus the inverse of the Jacobian mapping times the determinant
%    jac        :  Determinant of the Jacobian mapping 

ngv   = size(shapgeomvt,1);   % number of gauss points
npv   = size(shapgeomvt,2);
nd    = size(dgnodes,2);
ne    = size(dgnodes,3);

dshapvt   = reshape(permute(shapgeomvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
shapgeomvt    = shapgeomvt(:,:,1); 

% compute dg nodes at Gauss points: x = sum_i phi_i(xi) *x_i
pn = reshape(dgnodes,[npv nd*ne]);
pg = shapgeomvt*pn;
pg = reshape(pg,[ngv nd ne]);

% compute the Jacobian matrix at Gauss points: dx/dxi
Jg = dshapvt*pn;
Jg = permute(reshape(Jg,[ngv nd nd ne]),[2 3 1 4]);
Jg = reshape(Jg,[nd,nd,ngv*ne]);






