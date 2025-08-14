function [A, B1, B2, C1, C2, D1, D2, DL, DU] = pathcompute2(AE, epath, fpath, lpath, fintf, lintf, f2e, e2f, elcon)

[A, B1, B2, C1, C2, D1, D2, DL, DU] = pathsystem2(AE, epath, fpath, lpath, fintf, lintf, f2e, e2f, elcon);

ne = size(A,3);
for i = 1:ne
    A(:,:,i) = inv(A(:,:,i));    
    C1(:,:,i) = C1(:,:,i)*A(:,:,i);    
    C2(:,:,i) = C2(:,:,i)*A(:,:,i);    
end

m = size(DL, 1);
[npaths, nep] = size(epath);

q = size(A,1);
A = reshape(A,[q q npaths nep]);
B1 = reshape(B1,[q m npaths nep]);
B2 = reshape(B2,[q m npaths nep]);
C1 = reshape(C1,[m q npaths nep]);
C2 = reshape(C2,[m q npaths nep]);
D1 = reshape(D1,[m m npaths nep]);
D2 = reshape(D2,[m m npaths nep]);
DL = reshape(DL,[m m npaths nep]);
DU = reshape(DU,[m m npaths nep]);

for k = 1:nep
  for n = 1:npaths
    CB11 = C1(:,:,n,k)*B1(:,:,n,k);
    CB12 = C1(:,:,n,k)*B2(:,:,n,k);
    CB21 = C2(:,:,n,k)*B1(:,:,n,k);
    CB22 = C2(:,:,n,k)*B2(:,:,n,k);
    D1(:,:,n,k) = D1(:,:,n,k) - CB11;
    D2(:,:,n,k) = D2(:,:,n,k) - CB22;
    DU(:,:,n,k) = DU(:,:,n,k) - CB12;
    DL(:,:,n,k) = DL(:,:,n,k) - CB21;
    if k < nep        
      D1(:,:,n,k+1) = D1(:,:,n,k+1) - CB22; 
    end  
  end
end

% LU decomposition of tridagonal system
for k = 1:nep
  for n = 1:npaths
    D1(:,:,n,k) = inv(D1(:,:,n,k));
    DL(:,:,n,k) = DL(:,:,n,k)*D1(:,:,n,k);   
    if k < nep
      D1(:,:,n,k+1) = D1(:,:,n,k+1) - DL(:,:,n,k)*DU(:,:,n,k);
    else
      D2(:,:,n,k) = D2(:,:,n,k) - DL(:,:,n,k)*DU(:,:,n,k);
      D2(:,:,n,k) = inv(D2(:,:,n,k));
    end
  end  
end


