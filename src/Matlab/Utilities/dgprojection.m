function UDG1 = dgprojection(master1,mesh1,UDG,porder)

nc    = size(UDG,2);
nd    = master1.nd;
npv1  = master1.npv;
ngv1  = master1.ngv;
ne    = size(mesh1.dgnodes,3);

dshapvt1  = reshape(permute(master1.shapvt(:,:,2:nd+1),[1 3 2]),[ngv1*nd npv1]);
plocvl = masternodes(porder,nd,mesh1.elemtype,mesh1.nodetype);
shapt = mkshape(porder,plocvl,master1.gpvl,mesh1.elemtype);
shap  = squeeze(shapt(:,:,1));

UDG1 = zeros(npv1,nc,ne);
for i=1:ne
    dg1 = mesh1.dgnodes(:,:,i);        

    % compute the Jacobian matrix at Gauss points: dx/dxi    
    Jg1 = dshapvt1*dg1(:,1:nd);
    Jg1 = reshape(Jg1,[ngv1 nd nd]);        
    jac1 = volgeom(Jg1);        
            
    M1 = master1.shapvl(:,:,1)*diag(master1.gwvl.*jac1)*master1.shapvl(:,:,1)';
    C1 = master1.shapvl(:,:,1)*diag(master1.gwvl.*jac1)*shap';
    L1 = C1*UDG(:,:,i);
    UDG1(:,:,i) = M1\L1;        
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




