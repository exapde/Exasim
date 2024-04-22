function UH=inituhat(master,elcon,UDG,ncu)

if isstruct(master)
    perm  = master.perm(:,:,1);
else
    perm = master(:,:,1);
end

ne   = size(UDG,3);
npf  = size(perm,1);
nfe  = size(perm,2);
nsiz = max(elcon(:));

UH   = zeros(ncu*nsiz,1);
ndiv = zeros(nsiz,1);

for i=1:ne
    con = repmat((elcon(:,i)'-1)*ncu,ncu,1)+repmat((1:ncu)',1,nfe*npf);
    il = reshape(con,nfe*npf*ncu,1);
    fe = reshape(permute(UDG(reshape(perm,nfe*npf,1),1:ncu,i),[2,1]),ncu*npf*nfe,1);
    UH(il) = UH(il) + fe;
    ndiv(elcon(:,i)) = ndiv(elcon(:,i)) + 1;
end

ndiv = 1./ndiv;
UH = reshape(UH,ncu,nsiz);
UH = bsxfun(@times,UH,reshape(ndiv,[1 nsiz]));

