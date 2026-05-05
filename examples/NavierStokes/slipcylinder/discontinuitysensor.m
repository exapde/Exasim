function [sensor, s2] = discontinuitysensor(master, mesh, s)

nd = mesh.nd;
porder = 2;
if nd == 1
    np1 = porder;
    A = tensorproduct(master.xpe, mesh.porder);  
elseif nd==2
    if mesh.elemtype==0
    np1 = porder*(porder+1)/2;
    A = koornwinder(master.xpe, mesh.porder);  
    else
    np1 = porder*porder;
    A = tensorproduct(master.xpe, mesh.porder);
    end
else
    if mesh.elemtype==0
    np1 = porder*(porder+1)*(porder+2)/6;
    A = koornwinder(master.xpe, mesh.porder);  
    else
    np1 = porder*porder*porder;
    A = tensorproduct(master.xpe, mesh.porder);
    end
end

ne = size(mesh.dgnodes,3);
u = A\squeeze(s);
u((np1+1):end,:) = 0;
U = reshape(A*u,[size(mesh.dgnodes,1) 1 ne]); 
s2 = calerror(mesh,master,s,U,2);
 
sensor = 0*mesh.dgnodes(:,1,:);
for i = 1:ne
    sensor(:,1,i) = s2(i);
end

end

function err = calerror(mesh,master,UDG, U0, p)
    [npv, nc, ne] = size(UDG);
    ngv = master.ngv;
    nd = master.nd;
    shapvt = squeeze(master.shapvt(:,:,1));
    dshapvt = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);
    err = zeros(nc,ne);
    for i = 1:ne
        dg = mesh.dgnodes(:,:,i);
        % compute the Jacobian matrix at Gauss points: dx/dxi
        Jg = dshapvt*dg(:,1:nd);
        Jg = reshape(Jg,[ngv nd nd]); 
        jac = volgeom(Jg); 
        udgg = shapvt*UDG(:,:,i); 
        udge = shapvt*U0(:,:,i); 
        for j = 1:nc
            err(j,i) = (master.gwe.*jac)'*(abs(udgg(:,j)./udge(:,j) - 1).^p)/sum(master.gwe.*jac);     
        end  
    end
    err = sqrt(err);
end

function [jac] = volgeom(Jg)
    nd = size(Jg,2);
    switch nd
    case 1
    jac = Jg; 
    case 2
    jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);
    case 3
    jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+...
    Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+...
    Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);
    otherwise
    error('Dimension is not implemented');
    end
end

