function L = averagevector(master,mesh)

nd    = master.nd;
npv   = master.npv;
ngv   = master.ngv;
ne    = size(mesh.dgnodes,3);
dshapegt   = reshape(permute(master.shapegt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
L = zeros(npv,1,ne);
for i=1:ne
    dg  = mesh.dgnodes(:,:,i);
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapegt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    [jac,~] = volgeom(Jg);            
    L(:,1,i)  = master.shapegw(:,:,1)*(jac);
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

end




