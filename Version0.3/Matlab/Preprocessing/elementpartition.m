function dmd = elementpartition(t,elemtype,nproc,metis)
 
[nve,ne] = size(t);
dmd = cell(nproc,1);

if nproc==1
    i = 1;
    dmd{i}.elempart = 1:ne;
    dmd{i}.elempartpts = ne;
    dmd{i}.nbsd = [];
    dmd{i}.elem2cpu = [];
    dmd{i}.elemrecv = [];    
    dmd{i}.elemsend = [];
    dmd{i}.elemrecvpts = [];    
    dmd{i}.elemsendpts = [];
    return;
end

elem2cpu = partition(t',ne,nproc,metis);

if nve==2
    dim=1;    
else
    if elemtype==0 % tri/tet elements
        dim=nve-1;            
    elseif elemtype==1 % quad/hex elements
        dim=log2(nve);           
    end
end
face = getelemface(dim,elemtype);
[nvf,nfe] = size(face);

disp('run mkv2t...');  
[re,ce] = mkv2t(t,ne);

for i = 1:nproc
    disp(['preprocessing ' num2str(i)]); 
    intelem = find(elem2cpu == (i-1)); % elements in subdomain i    
    elem = node2elem(t(:,intelem),re,ce);    
    extelem = sort(setdiff(elem,intelem)); % exterior elements    
    
    % fix extelem
    n1 = length(intelem);
    n2 = length(extelem);    
    t1 = reshape(t(face,intelem),[nvf nfe*n1]);
    t2 = reshape(t(face,extelem),[nvf nfe n2]);    
    t1 = sort(t1,1);
    t2 = sort(t2,1);
    match = zeros(1,n2);
    for j=1:n2        
        for k=1:nfe
            if min(sum(abs(t2(:,k,j) - t1),1))==0
                match(j) = 1;
                break;
            end
        end                
    end
    extelem = extelem(match==1);
            
    elem = node2elem(t(:,extelem),re,ce); % all elements connected to exterior elements
    bndelem = intersect(elem,intelem);  % boundary elements in subdmain i
    outelem = sort(setdiff(elem,[intelem; extelem])); %  elements outside subdmain i            
    
    % fix outelem
    n1 = length([intelem; extelem]);
    n2 = length(outelem);    
    t1 = reshape(t(face,[intelem; extelem]),[nvf nfe*n1]);
    t2 = reshape(t(face,outelem),[nvf nfe n2]);    
    t1 = sort(t1,1);
    t2 = sort(t2,1);
    match = zeros(1,n2);
    for j=1:n2        
        for k=1:nfe
            if min(sum(abs(t2(:,k,j) - t1),1))==0
                match(j) = 1;
                break;
            end
        end                
    end    
    outelem = outelem(match==1);
        
    dmd{i}.elempart = [setdiff(intelem,bndelem); bndelem; extelem; outelem]; % partitions of elements
    dmd{i}.elempartpts = [length(intelem)-length(bndelem) length(bndelem) length(extelem) length(outelem)];    
    dmd{i}.elem2cpu = elem2cpu(dmd{i}.elempart);
    nelem = dmd{i}.elempartpts;

    recvelem = [extelem; outelem]; % elements received from neighbors
    ind = xiny(recvelem(:), dmd{i}.elempart);    
    dmd{i}.elemrecv = [dmd{i}.elem2cpu(ind)+1 (sum(nelem(1:2))+1:sum(nelem))' recvelem];        
    [~,ind] = sort(dmd{i}.elemrecv(:,1)); 
    dmd{i}.elemrecv = dmd{i}.elemrecv(ind,:);
    dmd{i}.nbsd = unique(dmd{i}.elemrecv(:,1)); % neighboring subdomains         
end

% store elements sent to neighboring subdomains to assemble the linear system
for k = 1:nproc    
    dmd{k}.elemsend = [];
end
for i = 1:nproc            
    for j = 1:length(dmd{i}.nbsd)
        % cpu k sends information to cpu i
        k = dmd{i}.nbsd(j);
        ii = dmd{i}.elemrecv(:,1)==k;
        tm = dmd{i}.elemrecv(ii,:);
        tm(:,1) = i;        
        tm(:,2) = xiny(tm(:,3),dmd{k}.elempart(:));
        dmd{k}.elemsend = [dmd{k}.elemsend; tm];        
    end    
end

for i = 1:nproc    
    for j = 1:length(dmd{i}.nbsd)
        dmd{i}.elemsendpts(j) = length(find(dmd{i}.elemsend(:,1)==dmd{i}.nbsd(j)));
        dmd{i}.elemrecvpts(j) = length(find(dmd{i}.elemrecv(:,1)==dmd{i}.nbsd(j)));
    end
    dmd{i}.elemsend = dmd{i}.elemsend(:,2);
    dmd{i}.elemrecv = dmd{i}.elemrecv(:,2);
end

% elempart = cell(nproc,1);
% elempartall = cell(nproc,1);
% for i = 1:nproc
%     nei = sum(dmd{i}.elempartpts(1:2));
%     elempart{i} = dmd{i}.elempart(1:nei);
%     elempartall{i} = dmd{i}.elempart;
% end
% save dmd.mat elempart elempartall;

