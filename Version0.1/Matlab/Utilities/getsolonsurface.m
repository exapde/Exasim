function [u,dgnodes] = getsolonsurface(mesh,UDG,ib)

perm = mesh.perm;
nf = mesh.nf;
nd = mesh.nd;
npf = size(perm,1);

in = find(mesh.f(:,end)==ib);
ns = length(in); 
dgnodes = zeros(npf,nd,ns);

nc = size(UDG,2);
u = zeros(npf,nc,ns);
for j = 1:ns
    i = in(j);
    fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
    kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
    i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element                            
    p = mesh.dgnodes(perm(:,i1),1:nd,fi(1));        
    dgnodes(:,:,j) = p;
    u(:,:,j) = UDG(perm(:,i1),:,fi(1));                                
end

% in = ismember(mesh.f(:,end),ib);
% ns = sum(in); 
% dgnodes = zeros(npf,nd,ns);
% 
% nc = size(UDG,2);
% u = zeros(npf,nc,ns);
% j = 1;
% for i = 1:nf
%     fi = mesh.f(i,end-1:end); % obtain two elements sharing the same face i  
%     if ismember(fi(2),ib) % face i is a boundary face
%         kf = mesh.t2f(fi(1),:);    % obtain neighboring faces 
%         i1 = find(kf(1,:)==i);  % obtain the index of face i in the 1st element                            
%         p = mesh.dgnodes(perm(:,i1),1:nd,fi(1));        
%         dgnodes(:,:,j) = p;
%         u(:,:,j) = UDG(perm(:,i1),:,fi(1));                
%         j = j + 1;                
%     end
% end

