function [err,rerr] = computeerror(UDG,mesh,master,func,time)

[npe,nc,ne] = size(UDG);
nge = master.nge;
nd  = master.dim;

shapvt    = squeeze(master.shapegt(:,:,1));
dshapvt   = reshape(permute(master.shapegt(:,:,2:nd+1),[1 3 2]),[nge*nd npe]);

err = zeros(nc,1);
errU = zeros(nc,1);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[nge nd nd]);        
    jac = volgeom(Jg);
                
    pg = shapvt*dg;    
    if nargin<5
        udge = feval(func,pg);
    else
        udge = feval(func,pg,time);
    end
    
    udgg = shapvt*UDG(:,:,i);    
    for j = 1:nc
        err(j) = err(j) + (master.gwe.*jac)'*(udgg(:,j)-udge(:,j)).^2;  
        errU(j) = errU(j) + (master.gwe.*jac)'*udge(:,j).^2;  
    end    
end
err  = sqrt(err);
rerr = err./sqrt(errU);

function [jac] = volgeom(Jg)

nd  = size(Jg,2);
switch nd
    case 1
        jac = Jg;        
    case 2
        jac = Jg(:,1,1).*Jg(:,2,2) - Jg(:,1,2).*Jg(:,2,1);        
    case 3
        jac = Jg(:,1,1).*Jg(:,2,2).*Jg(:,3,3) - Jg(:,1,1).*Jg(:,3,2).*Jg(:,2,3)+ ...
              Jg(:,2,1).*Jg(:,3,2).*Jg(:,1,3) - Jg(:,2,1).*Jg(:,1,2).*Jg(:,3,3)+ ...
              Jg(:,3,1).*Jg(:,1,2).*Jg(:,2,3) - Jg(:,3,1).*Jg(:,2,2).*Jg(:,1,3);                    
    otherwise
        error('Dimension is not implemented');
end




