function f = mkf(t,fb,nd)

[nve,ne] = size(t);

elemtype = 0;
if (nd==2) && (nve==4)
    elemtype=1;    
end
if (nd==3) && (nve==8)
    elemtype=1;    
end

face = getelemface(nd,elemtype);
[nvf,nfe] = size(face);
tf = reshape(t(face,:),[nvf nfe ne]);
f2e = mkf2e(t,elemtype,nd);

nf = size(f2e,2);
f = zeros(nvf+4,nf);
for i = 1:nf
    f(1:nvf,i) = tf(:,f2e(2,i),f2e(1,i));
end
f(nvf+1,:) = f2e(1,:);
f(nvf+2,:) = f2e(3,:);
f(nvf+3,:) = f2e(2,:);
f(nvf+4,:) = f2e(4,:);

ind = find(f(nvf+2,:)==0);
for i = 1:length(ind)
    e1 = f(nvf+1,ind(i));
    l1 = f(nvf+3,ind(i));
    f(nvf+4,ind(i)) = -fb(l1,e1);
end



