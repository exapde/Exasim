function facecon = faceconnectivity(p,t,f,f2t,elemtype,prdexpr,porder)

dim = size(p,1);
nprd = size(prdexpr,1);
nf = size(f2t,2);

[philocvl,philocfc,~,~,perm] = localbasis(porder,dim,elemtype);
philocvl = philocvl';
philocfc = philocfc';
face = getelemface(dim,elemtype);

npf = size(perm,1);
npe = size(philocvl,2);

facecon = zeros(npf,2,nf);
for i = 1:nf
    e1 = f2t(1,i);
    l1 = f2t(2,i);    
    e2 = f2t(3,i);
    l2 = f2t(4,i);    
            
    f1 = t(face(:,l1),e1);    
    pf = p(:,f1)*philocfc; % nodal points on face i    
    pf = round(pf,8);    
        
    t1 = t(:,e1);
    p1 = p(:,t1)*philocvl;
    p1 = p1(:,perm(:,l1));
    p1 = round(p1,8);        
    
    if f(l1,e1)<0 % periodic face        
        pbnd = -f(l1,e1);
        for k=1:nprd
            if (prdexpr{k,1}==pbnd)
                pf = prdexpr{k,2}(pf);
                p1 = prdexpr{k,2}(p1);
                break;
            elseif (prdexpr{k,3}==pbnd)
                pf = prdexpr{k,4}(pf);
                p1 = prdexpr{k,4}(p1);
                break;
            end
        end        
    end
    
    j1 = xiny(p1',pf');  
    facecon(:,1,i) = (e1-1)*npe + perm(j1,l1);      
    
    if e2>0 % face i is an interior face            
        t2 = t(:,e2);
        p2 = p(:,t2)*philocvl;      
        p2 = p2(:,perm(:,l2));
        p2 = round(p2,8);                                       
        if f(l2,e2)<0 % periodic face        
            pbnd = -f(l2,e2);                        
            if (prdexpr{k,1}==pbnd)
                p2 = prdexpr{k,2}(p2);                    
            elseif (prdexpr{k,3}==pbnd)
                p2 = prdexpr{k,4}(p2);                    
            end                            
        end        
        j2 = xiny(p2',pf');     
        facecon(:,2,i) = (e2-1)*npe + perm(j2,l2);       
    else
        facecon(:,2,i) = facecon(:,1,i);
    end    
end 




