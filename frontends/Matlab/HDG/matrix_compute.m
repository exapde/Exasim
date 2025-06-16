function [A1, A2, A3, C1, C2, C3, D] = matrix_compute(A1, A2, A3, B1, B2, B3, C1, C2, C3, D, idx1, idx2, idx3)

M = size(A1,1);
nb = size(A1, 3);
count1 = size(A1,4);
count2 = size(A2,4);
count3 = size(A3,4);

for i = 1:count1
  for n = 1:nb
    A1(:,:,n,i) = inv(A1(:,:,n,i));
  end
end
for i = 1:count2
  for n = 1:nb
    A2(:,:,n,i) = inv(A2(:,:,n,i));
  end
end
for i = 1:count3
  for n = 1:nb
    A3(:,:,n,i) = inv(A3(:,:,n,i));
  end
end

for i = 1:count1
  for n = 1:nb
    C1(:,:,n,i) = C1(:,:,n,i)*A1(:,:,n,i);
  end
end
for i = 1:count2
  for n = 1:nb
    C2(:,:,n,i) = C2(:,:,n,i)*A2(:,:,n,i);
  end
end
for i = 1:count3
  for n = 1:nb
    C3(:,:,n,i) = C3(:,:,n,i)*A3(:,:,n,i);
  end
end

N = size(idx1,1);
for i = 1:count1
  for j = 1:nb
    CAB = reshape(C1(:,:,j,i)*B1(:,:,j,i), [M N M N]);
    for m = 1:N
      for n = 1:N
        k = idx1(m, n, i);
        D(:,:,j,k) = D(:,:,j,k) - reshape(CAB(:,m,:,n),[M M]);
      end
    end
  end
end

N = size(idx2,1);
for i = 1:count2
  for j = 1:nb
    CAB = reshape(C2(:,:,j,i)*B2(:,:,j,i), [M N M N]);
    for m = 1:N
      for n = 1:N
        k = idx2(m, n, i);
        D(:,:,j,k) = D(:,:,j,k) - reshape(CAB(:,m,:,n), [M M]);
      end
    end
  end
end

N = size(idx3,1);
for i = 1:count3
  for j = 1:nb
    CAB = reshape(C3(:,:,j,i)*B3(:,:,j,i), [M N M N]);  
    for m = 1:N
      for n = 1:N
        k = idx3(m, n, i);
        D(:,:,j,k) = D(:,:,j,k) - CAB(:,m,:,n);
      end
    end
  end
end







