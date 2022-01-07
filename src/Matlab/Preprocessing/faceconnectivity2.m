function facecon = faceconnectivity2(t,f2t,dim,elemtype,porder)

nf = size(f2t,2);

[philocvl,~,~,plocfc,perm] = localbasis(porder,dim,elemtype);
permind = permindex(plocfc,dim,elemtype);
face = getelemface(dim,elemtype);

npf = size(perm,1);
npe = size(philocvl,1);

facecon = zeros(npf,2,nf);
if dim<=2
    for i = 1:nf
        e1 = f2t(1,i);
        l1 = f2t(2,i);    
        e2 = f2t(3,i);
        l2 = f2t(4,i);    

        facecon(:,1,i) = (e1-1)*npe + perm(:,l1);   
        if e2>0 % face i is an interior face                  
            facecon(:,2,i) = (e2-1)*npe + perm(permind,l2);                     
        else
            facecon(:,2,i) = facecon(:,1,i);
        end            
    end 
else
    for i = 1:nf
        e1 = f2t(1,i);
        l1 = f2t(2,i);    
        e2 = f2t(3,i);
        l2 = f2t(4,i);    

        facecon(:,1,i) = (e1-1)*npe + perm(:,l1);   
        if e2>0 % face i is an interior face           
            f1 = t(face(:,l1),e1);    
            f2 = t(face(:,l2),e2);                          
            if elemtype==0
                if (f1(1)==f2(1)) && (f1(2)==f2(3)) && (f1(3)==f2(2))
                    k = 1;
                elseif (f1(1)==f2(2)) && (f1(2)==f2(1)) && (f1(3)==f2(3))
                    k = 2;
                elseif (f1(1)==f2(3)) && (f1(2)==f2(2)) && (f1(3)==f2(1))    
                    k = 3;
                else
                    error("Mesh connectivity is wrong");
                end
            else                
                if (f1(1)==f2(1)) && (f1(2)==f2(4)) && (f1(3)==f2(3)) && (f1(4)==f2(2))
                    k = 1;
                elseif (f1(1)==f2(4)) && (f1(2)==f2(3)) && (f1(3)==f2(2)) && (f1(4)==f2(1))
                    k = 2;
                elseif (f1(1)==f2(3)) && (f1(2)==f2(2)) && (f1(3)==f2(1)) && (f1(4)==f2(4))
                    k = 3;
                elseif (f1(1)==f2(2)) && (f1(2)==f2(1)) && (f1(3)==f2(4)) && (f1(4)==f2(3))
                    k = 4;    
                else
                    error("Mesh connectivity is wrong");
                end
            end
            facecon(:,2,i) = (e2-1)*npe + perm(permind(:,k),l2);                     
        else
            facecon(:,2,i) = facecon(:,1,i);
        end    
    end     
end
