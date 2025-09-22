function q = qequationschur(MiC, MiE, u, uhat, s, fc)
% fc * M*q - C*u + E*uhat = M*s
% fc * q - Minv * C * u - Minv * E * uhat = s
% q = (1/fc) * (s + Minv * C * u - E * uhat) 

uhat = permute(uhat,[2 1 3]);

[npe, ncu, ne] = size(u);
nd = size(MiC,4);
q = reshape(s, [npe ncu nd ne]);
for i = 1:ne    
    for l = 1:nd            
      q(:,:,l,i) = q(:,:,l,i) + MiC(:,:,i,l)*u(:,:,i) - MiE(:,:,i,l)*uhat(:,:,i);     
    end    
end
q = (1/fc) * reshape(q, [npe ncu*nd ne]);
