function ucg = dg2cg2(udg, cgelcon, colent2elem, rowent2elem)

[npe,ncu,ne] = size(udg);

nent = length(rowent2elem)-1;
vcg = zeros(nent, ncu);
for i = 1:nent
    elem = colent2elem((rowent2elem(i)+1):rowent2elem(i+1));        
    tm = udg(:,:,elem);
    vcg(i,:) = mean(reshape(tm,[npe*length(elem) ncu]),1);
end

ucg = zeros(npe,ncu,ne);
for i = 1:ne
    ucg(:,:,i) = vcg(cgelcon(:,i),:);
end
