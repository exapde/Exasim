function [A, B1, B2, C1, C2, D1, D2, DL, DU] = pathsystem2(AE, epath, fpath, lpath, fintf, lintf, f2e, e2f, elcon)

n1 = size(fintf,1);
[nfe,ne] = size(e2f);
npf = numel(elcon)/(nfe*ne);
ncu = sqrt(numel(AE)/(ne*npf*npf*nfe*nfe));

elcon = reshape(elcon, [npf nfe ne]);
AE = reshape(AE, [ncu npf nfe ncu npf nfe ne]);

epath = epath(:);
ne = length(epath(:));
A = zeros(ncu, npf, n1, ncu, npf, n1, ne);
B1 = zeros(ncu, npf, n1, ncu, npf, 1, ne);
B2 = zeros(ncu, npf, n1, ncu, npf, 1, ne);
C1 = zeros(ncu, npf, 1, ncu, npf, n1, ne);
C2 = zeros(ncu, npf, 1, ncu, npf, n1, ne);
D1 = zeros(ncu, npf, 1, ncu, npf, 1, ne);
D2 = zeros(ncu, npf, 1, ncu, npf, 1, ne);
DL = zeros(ncu, npf, 1, ncu, npf, 1, ne);
DU = zeros(ncu, npf, 1, ncu, npf, 1, ne);

for i0 = 1:ne
  for i1 = 1:n1
    e = epath(i0);
    f = fintf(i1,i0);
    e1 = f2e(1,f);
    l1 = f2e(2,f);  
    e2 = f2e(3,f);
    l2 = f2e(4,f);            
    j1 = elcon(:,l1,e1) - npf*(f-1);

    % diagonal blocks
    if e2>0
      j2 = elcon(:,l2,e2) - npf*(f-1);  
      A(:,:,i1,:,:,i1,i0) = AE(:,j1,l1,:,j1,l1,e1) + AE(:,j2,l2,:,j2,l2,e2);          
    else    
      A(:,:,i1,:,:,i1,i0) = AE(:,j1,l1,:,j1,l1,e1);
    end  
       
    l1 = lintf(i1,i0);
    j1 = elcon(:,l1,e) - npf*(f-1);     
    
    % off-sdiagonal blocks
    for i2 = 1:n1
      if i1 ~= i2
        l2 = lintf(i2,i0);
        f2 = fintf(i2,i0);
        j2 = elcon(:,l2,e) - npf*(f2-1);     
        A(:,:,i1,:,:,i2,i0) = AE(:,j1,l1,:,j2,l2,e);
      end
    end
            
    % form B
    m = 1;    
    f2 = fpath(m,i0);
    l2 = lpath(m,i0);
    j2 = elcon(:,l2,e) - npf*(f2-1);           
    B1(:,:,i1,:,:,1,i0) = AE(:,j1,l1,:,j2,l2,e);
    C1(:,:,1,:,:,i1,i0) = AE(:,j2,l2,:,j1,l1,e);    
    
    m = 2;    
    f2 = fpath(m,i0);
    l2 = lpath(m,i0);
    j2 = elcon(:,l2,e) - npf*(f2-1);           
    B2(:,:,i1,:,:,1,i0) = AE(:,j1,l1,:,j2,l2,e);
    C2(:,:,1,:,:,i1,i0) = AE(:,j2,l2,:,j1,l1,e);    
        
  end
end

for i0 = 1:ne
  e = epath(i0);
  
  f1 = fpath(1,i0);
  l1 = lpath(1,i0);
  j1 = elcon(:,l1,e) - npf*(f1-1);                 

  f2 = fpath(2,i0);
  l2 = lpath(2,i0);
  j2 = elcon(:,l2,e) - npf*(f2-1);                     
  
  D1(:,:,1,:,:,1,i0) = AE(:,j1,l1,:,j1,l1,e);
  DU(:,:,1,:,:,1,i0) = AE(:,j1,l1,:,j2,l2,e);  
  DL(:,:,1,:,:,1,i0) = AE(:,j2,l2,:,j1,l1,e);     
  D2(:,:,1,:,:,1,i0) = AE(:,j2,l2,:,j2,l2,e);
  
  if e == f2e(1,f1)
    e3 = f2e(3,f1);
    l3 = f2e(4,f1);
  elseif e == f2e(3,f1)    
    e3 = f2e(1,f1);
    l3 = f2e(2,f1);
  end
      
  if e3>0
    j3 = elcon(:,l3,e3) - npf*(f1-1);      
    D1(:,:,1,:,:,1,i0) = D1(:,:,1,:,:,1,i0) + AE(:,j3,l3,:,j3,l3,e3);      
  end  
  
  if e == f2e(1,f2)
    e3 = f2e(3,f2);
    l3 = f2e(4,f2);
  elseif e == f2e(3,f2)    
    e3 = f2e(1,f2);
    l3 = f2e(2,f2);
  end
    
  if e3>0
    j3 = elcon(:,l3,e3) - npf*(f2-1);      
    D2(:,:,1,:,:,1,i0) = D2(:,:,1,:,:,1,i0) + AE(:,j3,l3,:,j3,l3,e3);      
  end     
end
  
A = reshape(A, ncu*npf*n1, ncu*npf*n1, ne);
B1 = reshape(B1, ncu*npf*n1, ncu*npf, ne);
B2 = reshape(B2, ncu*npf*n1, ncu*npf, ne);
C1 = reshape(C1, ncu*npf, ncu*npf*n1, ne);
C2 = reshape(C2, ncu*npf, ncu*npf*n1, ne);
D1 = reshape(D1, ncu*npf*1, ncu*npf*1, ne);
D2 = reshape(D2, ncu*npf*1, ncu*npf*1, ne);
DL = reshape(DL, ncu*npf*1, ncu*npf*1, ne);
DU = reshape(DU, ncu*npf*1, ncu*npf*1, ne);

