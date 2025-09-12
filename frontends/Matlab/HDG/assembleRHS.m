function FE = assembleRHS(FE, elcon)

[ndf, ne] = size(elcon);
ncu = size(FE,1) / ndf;

il = zeros(ndf*ncu,ndf*ncu,ne);
for i=1:ne
    con = repmat((elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,ndf);
    con = reshape(con,ndf*ncu,1);    
    il(:,:,i) = repmat(con ,1,ndf*ncu);
end

FE = sparse(reshape(il(:,1,:),ndf*ncu*ne,1),ones(ndf*ncu*ne,1),reshape(FE,ndf*ncu*ne,1));                                                                      
