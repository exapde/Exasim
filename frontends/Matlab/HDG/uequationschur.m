function [AE, FE, DUDG, DUDG_DUH, D, F, K, H] = uequationschur(Ru, Rh, B, D, F, G, K, H, MinvC, MinvE)

[npe, ne, ncu] = size(Ru);
ndf = size(Rh,1);

if isempty(B) == 0
  nd = size(MinvC,1) / npe;
end

DUDG = zeros(npe*ncu, ne);
DUDG_DUH = zeros(npe*ncu, ncu*ndf, ne);
FE = zeros(ncu*ndf, ne);
AE = zeros(ncu*ndf, ncu*ndf, ne);

Ru = reshape(permute(Ru, [1 3 2]), [npe*ncu ne]); % Ru = npe * ne * ncu
Rh = reshape(permute(Rh, [3 1 2]), [ncu*ndf ne]); % Rh = ndf * ne * ncu     
D = reshape(permute(D, [1 4 2 5 3]), [npe*ncu npe*ncu ne]); % D = npe * npe * ne * ncu * ncu
F = reshape(permute(F, [1 4 5 2 3]), [npe*ncu ncu*ndf ne]); % F = npe * ndf * ne * ncu * ncu
K = reshape(permute(K, [4 1 2 5 3]), [ncu*ndf npe*ncu ne]); % K = ndf * npe * ne * ncu * ncu 
H = reshape(permute(H, [4 1 5 2 3]), [ncu*ndf ncu*ndf ne]); % H = ndf * ndf * ne * ncu * ncu 

if (isempty(B) == 0)                   
  % B = npe * npe * ne * ncu * ncu * nd
  B = reshape(permute(reshape(B,[npe npe ne ncu ncu nd]), [1 4 5 2 6 3]), [npe*ncu*ncu npe*nd ne]);       

  % G = ndf * npe * ne * ncu * ncu * nd  
  G = reshape(permute(reshape(G, [ndf npe ne ncu ncu nd]), [4 1 5 2 6 3]), [ncu*ndf*ncu npe*nd ne]);       

  for i = 1:ne        
     % MinvC = (npe * nd) * npe * ne 
     BMiC = reshape(B(:,:,i)*MinvC(:,:,i), [npe ncu ncu npe]);
     D(:,:,i) = D(:,:,i) + reshape(permute(BMiC, [1 2 4 3]),[npe*ncu npe*ncu]);

     % MinvE = (npe * nd) * ndf * ne
     BMiE = reshape(B(:,:,i)*MinvE(:,:,i), [npe*ncu ncu*ndf]);
     F(:,:,i) = F(:,:,i) - BMiE;

     GMiC = reshape(G(:,:,i)*MinvC(:,:,i), [ncu ndf ncu npe]);
     K(:,:,i) = K(:,:,i) + reshape(permute(GMiC, [1 2 4 3]), [ncu*ndf npe*ncu]);

     GMiE = reshape(G(:,:,i)*MinvE(:,:,i), [ncu*ndf ncu*ndf]);
     H(:,:,i) = H(:,:,i) - GMiE;
  end
end

for i = 1:ne          
  DUDG(:,i) = D(:,:,i)\Ru(:,i);          
  DUDG_DUH(:,:,i) = -D(:,:,i)\F(:,:,i);            
  FE(:,i) = Rh(:,i) - K(:,:,i)*DUDG(:,i);
  AE(:,:,i) = H(:,:,i) + K(:,:,i)*DUDG_DUH(:,:,i);
end

