function qdg = gradu(shapen, dgnodes, udg)
% shapen: npe*npe*nd, nodal*shape*dimension

% udg: npe*ne*ncu;  shapen: npe*npe*nd; Xx: npe*ne*nd*nd 
% shapen*udg: npe*ne*ncu*nd; shapen*udg*Xx: npe*ne*ncu*nd
    
[npe,ncu,ne] = size(udg);
nd = size(shapen,3);

udg = permute(udg,[1 3 2]);
Xx = volgeom(shapen, dgnodes);
tmp = zeros(npe,ne*ncu,nd);
for i = 1:nd
    tmp(:,:,i) = shapen(:,:,i)*reshape(udg,[npe ne*ncu]);
end
tmp = reshape(tmp,[npe*ne,ncu,nd]);

qdg = zeros(npe*ne,ncu,nd);
for j = 1:ncu
    for k = 1:nd
        for i = 1:nd        
            qdg(:,j,k) = qdg(:,j,k) + tmp(:,j,i).*Xx(:,k,i);
        end
    end
end
qdg = reshape(qdg,[npe ne ncu*nd]);
qdg = permute(qdg, [1 3 2]);

function Xx = volgeom(shapen,dgnodes)

[npe,nd,ne] = size(dgnodes);
dgnodes = permute(dgnodes, [1 3 2]);

dshapvt   = reshape(permute(shapen,[1 3 2]),[npe*nd npe]);

% compute the Jacobian matrix at Gauss points: dx/dxi
Jg = dshapvt*reshape(dgnodes,[npe ne*nd]);
Jg = permute(reshape(Jg,[npe nd ne nd]),[1 3 2 4]);
Jg = reshape(Jg,[npe*ne,nd,nd]);

% jac: the determinant of the Jacobian matrix at Gauss points
% Xx = inv(Jg)*jac: the inverse of the Jacobian matrix times 
%                   the determinant of the Jacobian matrix 
switch nd
    case 1
        jac = Jg;
        Xx = ones(npe*ne,1);
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
        Xx(:,1,1) = Jg(:,2,2);  % dxi/dx
        Xx(:,2,1) = -Jg(:,2,1); % dxi/dy
        Xx(:,1,2) = -Jg(:,1,2); % deta/dx
        Xx(:,2,2) = Jg(:,1,1);  % deta/dy
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);            
        Xx(:,1,1) = Jg(:,2,2).*Jg(:,3,3) - Jg(:,2,3).*Jg(:,3,2);
        Xx(:,2,1) = Jg(:,2,3).*Jg(:,3,1) - Jg(:,2,1).*Jg(:,3,3);
        Xx(:,3,1) = Jg(:,2,1).*Jg(:,3,2) - Jg(:,2,2).*Jg(:,3,1);
        Xx(:,1,2) = Jg(:,1,3).*Jg(:,3,2) - Jg(:,1,2).*Jg(:,3,3);
        Xx(:,2,2) = Jg(:,1,1).*Jg(:,3,3) - Jg(:,1,3).*Jg(:,3,1);
        Xx(:,3,2) = Jg(:,1,2).*Jg(:,3,1) - Jg(:,1,1).*Jg(:,3,2);
        Xx(:,1,3) = Jg(:,1,2).*Jg(:,2,3) - Jg(:,1,3).*Jg(:,2,2);
        Xx(:,2,3) = Jg(:,1,3).*Jg(:,2,1) - Jg(:,1,1).*Jg(:,2,3);
        Xx(:,3,3) = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    otherwise
        error('Dimension is not implemented');
end

Xx = Xx./jac;

