function f2e = mkf2e(elcon, f, t2f, perm)

np = max(perm(:));
nf = size(f,1);
ne = size(t2f,1);
[npf,nfe] = size(perm);
elcon = reshape(elcon,[npf nfe ne]);

f2e = zeros(npf,nf,2);
for i = 1:nf
    fi = f(i,end-1:end); % obtain two elements sharing the same face i      
    if fi(2)>0           % face i is an interior face                
        kf = t2f(fi,:);         % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element
        i2 = kf(2,:)==i;  % obtain the index of face i in the 2nd element                                            
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;        
        f2e(:,i,1) =  (fi(1)-1)*np+perm(j1,i1);
        f2e(:,i,2) =  (fi(2)-1)*np+perm(j2,i2);        
    else % face i is a boundary face
        kf = t2f(fi(1),:); % obtain neighboring faces 
        i1 = kf(1,:)==i;  % obtain the index of face i in the 1st element                 
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;       
        f2e(:,i,1) =  (fi(1)-1)*np+perm(j1,i1);        
        f2e(:,i,2) = f2e(:,i,1);
    end        
end 
f2e = reshape(f2e,[npf*nf 2])'; 

