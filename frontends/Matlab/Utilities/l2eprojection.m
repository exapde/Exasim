function UDG = l2eprojection(mesh,master,func,param,time,nc)

ne = mesh.ne;
npv = master.npv;
ngv = master.ngv;
nd  = master.nd;

shapvt    = squeeze(master.shapvt(:,:,1));
dshapvt   = reshape(permute(master.shapvt(:,:,2:nd+1),[1 3 2]),[ngv*nd npv]);

UDG = zeros(npv,nc,ne);
for i = 1:ne
    dg = mesh.dgnodes(:,:,i);
    
    % compute the Jacobian matrix at Gauss points: dx/dxi
    Jg = dshapvt*dg(:,1:nd);
    Jg = reshape(Jg,[ngv nd nd]);        
    jac = volgeom(Jg);
                
    pg = shapvt*dg;           
    fg = feval(func,pg,param,time);    
    
    M = reshape(master.shapvgdotshapvl(:,:,1)*jac,[npv npv]);        
    for j = 1:nc
        F = master.shapvg(:,:,1)*(jac.*fg(:,j));
        UDG(:,j,i) = M\F;        
    end    
end

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




