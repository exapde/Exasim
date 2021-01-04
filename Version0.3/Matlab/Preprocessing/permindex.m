function ind = permindex(plocfc,dim,elemtype)

npf = size(plocfc,1);

if dim==1
    ind = 1;
elseif dim==2
    ind = (npf:-1:1)';
elseif dim==3
    if elemtype==0
        ind = zeros(npf,3);
        
        % [1 3 2]
        plocfc2 = plocfc; 
        plocfc2(:,1) = plocfc(:,2);
        plocfc2(:,2) = plocfc(:,1);
        ind(:,1) = xiny(round(plocfc,8),round(plocfc2,8));
    
        % [2 1 3]
        plocfc2 = plocfc; 
        plocfc2(:,1) = 1-plocfc(:,1)-plocfc(:,2);
        ind(:,2) = xiny(round(plocfc,8),round(plocfc2,8));            

        % [3 2 1]
        plocfc2 = plocfc;         
        plocfc2(:,2) = 1-plocfc(:,1)-plocfc(:,2);
        ind(:,3) = xiny(round(plocfc,8),round(plocfc2,8));                            
    else
        ind = zeros(npf,4);

        % [1 4 3 2]
        plocfc2 = plocfc; 
        plocfc2(:,1) = plocfc(:,2);
        plocfc2(:,2) = plocfc(:,1);
        ind(:,1) = xiny(round(plocfc,8),round(plocfc2,8));

        % [2 1 4 3]
        plocfc2 = plocfc; 
        plocfc2(:,2) = 1-plocfc(:,2);
        ind(:,2) = xiny(round(plocfc,8),round(plocfc2,8));            

        % [3 2 1 4]
        plocfc2 = plocfc; 
        plocfc2(:,1) = 1-plocfc(:,2);
        plocfc2(:,2) = 1-plocfc(:,1);
        ind(:,3) = xiny(round(plocfc,8),round(plocfc2,8));
        
        % [4 3 2 1]
        plocfc2 = plocfc;         
        plocfc2(:,1) = 1-plocfc(:,1);
        ind(:,4) = xiny(round(plocfc,8),round(plocfc2,8));                                    
    end
end
