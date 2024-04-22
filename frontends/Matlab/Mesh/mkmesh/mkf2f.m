function [f2f,f2t1,f2t2] = mkf2f(f, t2f)
%  
%   f     :  Face to element connectivity
%   t2f   :  Element to face connectivity
%
%   f2f   :  Face to face connectivity

nf  = size(f,1);    % number of faces
nfe = size(t2f,2);  % number of faces per element
nbf = 2*nfe-1;      % number of neighboring faces
f2f = zeros(nbf,nf);
f2t1 = zeros(nf,nfe+1);
f2t2 = zeros(nf,nfe+1);
for i = 1:nf
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    f2t1(i,1) = fi(1);
    f2t2(i,1) = fi(2);
    if fi(2)>0           % face i is an interior face
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element
        i2 = find(kf(2,:)==i);  % obtain the index of face i in the 2nd element            
%         f2t1(i,2:end) = [i1 setdiff(1:nfe,i1)];
%         f2t2(i,2:end) = [i2 setdiff(1:nfe,i2)];
        
        % the first block        
        k = 1;
        f2f(k,i) = i;        
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;                 
                j = kf(1,is);
                f2f(k,i) = j;
            end
        end
        
        for is=1:nfe % loop over faces of the 2nd element
            if is ~= i2                                                
                k = k + 1;                 
                j = kf(2,is);
                f2f(k,i) = j;
            end
        end        
    else % face i is a boundary face
        kf = t2f(fi(1),:);      % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element                 
%        f2t1(i,2:end) = [i1 setdiff(1:nfe,i1)];
        
        % the first block        
        k = 1;
        f2f(k,i) = i;        
        
        for is=1:nfe % loop over each faces of the 1st element
            if is ~= i1  
                k = k + 1;          
                j = kf(1,is);
                f2f(k,i) = j; 
            end
        end        
    end    
end

