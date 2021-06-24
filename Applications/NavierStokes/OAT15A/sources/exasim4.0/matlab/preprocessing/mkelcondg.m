function [elcon1, elcon2, udof] = mkelcondg(fblks,facecon,t2f,perm,npe)

elcon1 = [];
elcon2 = [];
n = size(fblks,2);    
udof = zeros(n+1,2);
for j=1:n
    f1 = fblks(1,j);
    f2 = fblks(2,j);
    elcon1 = [elcon1; mkelem2face(facecon,t2f,perm,npe,f1:f2,1)];
    elcon2 = [elcon2; mkelem2face(facecon,t2f,perm,npe,f1:f2,2)];    
    udof(j+1,1) = size(elcon1,1);
    udof(j+1,2) = size(elcon2,1);
end

function elcon = mkelem2face(facecon,t2f,perm,npe,facelist,opts)

nf = length(facelist);
npf = size(facecon,1);

edof = facecon(:,opts,facelist); % Global degrees of freedom of UDG on faces
edof = unique(edof); % remove duplications
n = length(edof);    % number of global degrees of freedom of UDG on faces
elcon = zeros(n,5);
elcon(:,1) = edof;

check = 1;
if check==1
    gdof = facecon(:,opts,facelist); gdof = gdof(:); % Global degrees of freedom of UDG on the minus-side faces
end

for i = 1:n % for each global degrees of freedom of UDG        
    elem = ceil(edof(i)/npe);       % element contains edof(i)
    eldof = edof(i) - (elem-1)*npe; % local index of edof(i)
    [~,lface] = find(perm==eldof);  % local faces connected to edof(i)             
    ind1 = [];
    for j = 1:length(lface) % for each local face
        fj = t2f(elem,lface(j));    % global face fj connected to edof(i)
        in = find(facelist==fj, 1); % index of global face fj in facelist
        if isempty(in)==0 % if global face fj is in facelist            
            ldof = find(facecon(:,opts,fj)==edof(i)); % find the 
            ind1 = [ind1 (in-1)*npf+ldof];
        end
    end
    elcon(i,2) = length(ind1);
    elcon(i,3:2+length(ind1)) = ind1;    
    
    if check
        ind = find(gdof==edof(i)); % find the corresponding local DOFs of UH on the minus-side faces        
        if max(abs(sort(ind1(:))-sort(ind(:))))~=0            
            error('something wrong');
        end
    end
end

% check 
ldof = elcon(:,3:end);
ldof = unique(ldof(:));
if ldof(1)==0
    ldof = ldof(2:end);
end
e = ldof-(1:npf*nf)';
if max(abs(e))>0
    error('Something wrong');
end
