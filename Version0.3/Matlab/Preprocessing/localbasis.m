function [philocvl,philocfc,plocvl,plocfc,perm] = localbasis(porder,dim,elemtype) 

[plocvl,~,plocfc,~,perm] = masternodes(porder,dim,elemtype);

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

