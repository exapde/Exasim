function mesh = fixdgmesh(mesh)

[mesh,nfixedfaces] = fixdgnodes(mesh);
str = ['number of fixed faces: ' num2str(nfixedfaces)]; disp(str);
while (nfixedfaces>0)
    [mesh,nfixedfaces] = fixdgnodes(mesh);
    str = ['number of fixed faces: ' num2str(nfixedfaces)]; disp(str);
end

function [mesh,nfixedfaces] = fixdgnodes(mesh)

perm = mesh.perm;
ne = mesh.ne;
nf = mesh.nf;
nd = mesh.nd;
[npf,nfe] = size(perm);
elcon = reshape(mesh.elcon,[npf nfe ne]);

nfixedfaces = 0;
for i = 1:nf
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
    if fi(2)>0        
        kf = mesh.t2f(fi(1),:);  % obtain neighboring faces 
        i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element        
        j1 = elcon(:,i1,fi(1)) - (i-1)*npf;                
        kf = mesh.t2f(fi(2),:);    % obtain neighboring faces 
        i2 = find(kf(1,:)==i);  % obtain the index of face i in the second element        
        j2 = elcon(:,i2,fi(2)) - (i-1)*npf;                            
        dg1 = mesh.dgnodes(perm(j1,i1),1:nd,fi(1));
        dg2 = mesh.dgnodes(perm(j2,i2),1:nd,fi(2));        
        dgm = (dg1-dg2);              
        if max(abs(dgm(:)))>1e-10
            nfixedfaces = nfixedfaces+1;
            %str = ['fix face ' num2str(i)]; disp(str);
            dgm = (dg1+dg2)/2;        
            mesh.dgnodes(perm(j1,i1),1:nd,fi(1)) = dgm;
            mesh.dgnodes(perm(j2,i2),1:nd,fi(2)) = dgm;      
        end
    end    
end


