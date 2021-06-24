function [philocvl,philocfc,plocvl,plocfc,perm] = localbasis(porder,dim,elemtype) 

[plocvl,~,plocfc,~,~,~,perm] = mkmasternodes(porder,dim,elemtype,0);

if dim==2 && elemtype==0      % tri
    nfe = length(plocfc);
    for i = 1:nfe        
        xi  = plocfc{i}(:,1);        
        philocfc{i}(:,1) = 1 - xi;
        philocfc{i}(:,2) = xi;
    end
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = 1 - xi - eta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;    
elseif dim==2 && elemtype==1  % quad
    nfe = length(plocfc);
    for i = 1:nfe        
        xi  = plocfc{i}(:,1);        
        philocfc{i}(:,1) = 1 - xi;
        philocfc{i}(:,2) = xi;
    end
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = (1-xi).*(1-eta);
    philocvl(:,2) = xi.*(1-eta);
    philocvl(:,3) = xi.*eta;
    philocvl(:,4) = (1-xi).*eta;
elseif dim==3 && elemtype==0  % tet
    nfe = length(plocfc);
    for i = 1:nfe        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);    
        philocfc{i}(:,1) = 1 - xi - eta;
        philocfc{i}(:,2) = xi;
        philocfc{i}(:,3) = eta;
    end
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = 1 - xi - eta - zeta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
    philocvl(:,4) = zeta;
elseif dim==3 && elemtype==1   % hex
    nfe = length(plocfc);
    for i = 1:nfe        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);
        philocfc{i}(:,1) = (1-xi).*(1-eta);
        philocfc{i}(:,2) = xi.*(1-eta);
        philocfc{i}(:,3) = xi.*eta;
        philocfc{i}(:,4) = (1-xi).*eta;
    end
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
elseif dim==3 && elemtype==2   % prism
    nfe = length(plocfc);
    for i = 1:2        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);    
        philocfc{i}(:,1) = 1 - xi - eta;
        philocfc{i}(:,2) = xi;
        philocfc{i}(:,3) = eta;
    end
    for i = 3:nfe        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);
        philocfc{i}(:,1) = (1-xi).*(1-eta);
        philocfc{i}(:,2) = xi.*(1-eta);
        philocfc{i}(:,3) = xi.*eta;
        philocfc{i}(:,4) = (1-xi).*eta;
    end    
    
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);    
    philocvl(:,1) = (1-xi-eta).*(1-zeta);
    philocvl(:,2) = xi.*(1-zeta);
    philocvl(:,3) = eta.*(1-zeta);      
    philocvl(:,4) = (1-xi-eta).*(zeta);
    philocvl(:,5) = xi.*(zeta);
    philocvl(:,6) = eta.*(zeta);    
elseif dim==3 && elemtype==3   % pyramid
    nfe = length(plocfc);
    for i = 1:1        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);
        philocfc{i}(:,1) = (1-xi).*(1-eta);
        philocfc{i}(:,2) = xi.*(1-eta);
        philocfc{i}(:,3) = xi.*eta;
        philocfc{i}(:,4) = (1-xi).*eta;
    end        
    for i = 2:nfe        
        xi  = plocfc{i}(:,1);
        eta = plocfc{i}(:,2);    
        philocfc{i}(:,1) = 1 - xi - eta;
        philocfc{i}(:,2) = xi;
        philocfc{i}(:,3) = eta;
    end
    
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);    
    philocvl(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocvl(:,2) = xi.*(1-eta).*(1-zeta);
    philocvl(:,3) = xi.*eta.*(1-zeta);
    philocvl(:,4) = (1-xi).*eta.*(1-zeta);    
    philocvl(:,5) = zeta;    
end

