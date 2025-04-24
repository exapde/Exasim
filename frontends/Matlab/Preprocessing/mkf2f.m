function [f2f, f2l] = mkf2f(f2e, e2f)
%  
%   f2e   :  Face to element connectivity
%   e2f   :  Element to face connectivity
%
%   f2f   :  Face to face connectivity

nf  = size(f2e,2);  % number of faces
nfe = size(e2f,1);  % number of faces per element
nbf = 2*(nfe-2);    % number of neighboring faces
f2f = zeros(nbf,nf);
f2l = zeros(nbf,nf);
for i = 1:nf  % loop over each face
    e1 = f2e(1,i);
    l1 = f2e(2,i);
    e2 = f2e(3,i);
    l2 = f2e(4,i);    
    
    k = 0;       
    for l=1:nfe % loop over each faces of the 1st element
        if l ~= l1  
            k = k + 1;                 
            j = e2f(l,e1);
            f2f(k,i) = j;
            f2l(k,i) = l;
        end
    end
    
    if e2>0           % face i is an interior face                
        for l=1:nfe % loop over faces of the 2nd element
            if l ~= l2                                                
                k = k + 1;                 
                j = e2f(l,e2);
                f2f(k,i) = j;
                f2l(k,i) = l;
            end
        end        
    end    
end


