function [r1, r2, r3, r4] = vector_compute(C1, C2, C3, rhs, idr1, idr2, idr3, nb)

%rhs = faceextract(b, face);

M = size(rhs, 1);
count1 = size(idr1, 2);
count2 = size(idr2, 2);
count3 = size(idr3, 2);

r1 = rhs(:,:,1:count1);
r2 = rhs(:,:,(count1+1):(count1+2*count2));
r3 = rhs(:,:,(count1+2*count2+1):(count1+2*count2+3*count3));
r4 = rhs(:,:,(count1+2*count2+3*count3+1):end);

r2 = reshape(permute(reshape(r2, [M nb 2 count2]),[1 3 2 4]), [M*2 nb count2]);
r3 = reshape(permute(reshape(r3, [M nb 3 count3]),[1 3 2 4]), [M*3 nb count3]);

if numel(C1)>0
  nfe = size(C1,1)/M + 1;
elseif numel(C2)>0
  nfe = size(C2,1)/M + 2;
elseif numel(C3)>0
  nfe = size(C3,1)/M + 3;  
end

for i = 1:count1
  for j = 1:nb
    ut = reshape(C1(:,:,j,i)*r1(:,j,i), [M nfe-1]);
    for k = 1:(nfe-1)
      m = idr1(k,i);
      r4(:,j,m) = r4(:,j,m) - ut(:,k);
    end
  end
end

for i = 1:count2
  for j = 1:nb
    ut = reshape(C2(:,:,j,i)*r2(:,j,i), [M nfe-2]);
    for k = 1:(nfe-2)
      m = idr2(k,i);
      r4(:,j,m) = r4(:,j,m) - ut(:,k);
    end
  end
end

for i = 1:count3
  for j = 1:nb
    ut = reshape(C3(:,:,j,i)*r3(:,j,i), [M nfe-3]);
    for k = 1:(nfe-3)
      m = idr3(k,i);
      r4(:,j,m) = r4(:,j,m) - ut(:,k);
    end
  end
end


