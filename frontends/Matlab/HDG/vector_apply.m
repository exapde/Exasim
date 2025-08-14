function uh = vector_apply(A1, A2, A3, B1, B2, B3, r1, r2, r3, u4, face, idr1, idr2, idr3)

count1 = size(idr1, 2);
count2 = size(idr2, 2);
count3 = size(idr3, 2);

[nb, nf] = size(face);
M = size(u4, 1);
uh = zeros(M, nb, nf);

if numel(B1)>0
  nfe = size(B1,2)/M + 1;
elseif numel(B2)>0
  nfe = size(B2,2)/M + 2;
elseif numel(B3)>0
  nfe = size(B3,2)/M + 3;  
end

ut = zeros(M, nfe-1);
for i = 1:count1
  for j = 1:nb   
    for k = 1:(nfe-1)
      m = idr1(k,i);
      ut(:,k) = u4(:,j,m);
    end    
    uh(:,j,i) = A1(:,:,j,i)*(r1(:,j,i) - B1(:,:,j,i)*ut(:));
  end
end

ut = zeros(M, nfe-2);
for i = 1:count2
  for j = 1:nb    
    for k = 1:(nfe-2)
      m = idr2(k,i);
      ut(:,k) = u4(:,j,m);            
    end
    u2 = A2(:,:,j,i)*(r2(:,j,i) - B2(:,:,j,i)*ut(:));
    q = count1 + 2*(i-1);
    uh(:,j,q + 1) = u2(1:M);
    uh(:,j,q + 2) = u2(M+1:2*M);
  end
end

ut = zeros(M, nfe-3);
for i = 1:count3
  for j = 1:nb
    for k = 1:(nfe-3)
      m = idr3(k,i);
      ut(:,k) = u4(:,j,m);            
    end
    u3 = A3(:,:,j,i)*(r3(:,j,i) - B3(:,:,j,i)*ut(:));
    q = count1 + count2*2 + 3*(i-1);
    uh(:,j,q + 1) = u3(1:M);
    uh(:,j,q + 2) = u3(M+1:2*M);
    uh(:,j,q + 3) = u3(2*M+1:3*M);
  end
end

q = count1 + count2*2 + count3*3 + 1;
uh(:,:,q:nf) = u4;

