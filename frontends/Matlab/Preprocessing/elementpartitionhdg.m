function dmd = elementpartitionhdg(t,t2t,nproc,metis)
 
[~,ne] = size(t);
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

for i = 1:nproc
    disp(['element partition ' num2str(i)]); 
    intelem = find(elem2cpu == (i-1)); % elements in subdomain i    
    elem = neighboringelements(t2t, intelem); % all elements connected to elements in subdomain i         
    extelem = sort(setdiff(elem,intelem)); % exterior elements    
                                
    elem = neighboringelements(t2t, extelem); % all elements connected to exterior elements
    bndelem = intersect(elem,intelem);  % interface elements in subdmain i
    
    dmd{i}.elempart = [setdiff(intelem,bndelem); bndelem; extelem]; % partitions of elements
    dmd{i}.elempartpts = [length(intelem)-length(bndelem) length(bndelem) length(extelem)];    
    dmd{i}.elem2cpu = elem2cpu(dmd{i}.elempart);
    nelem = dmd{i}.elempartpts;

    recvelem = extelem; % elements received from neighbors
    ind = xiny(recvelem(:), dmd{i}.elempart);    
    dmd{i}.elemrecv = [dmd{i}.elem2cpu(ind)+1 (sum(nelem(1:2))+1:sum(nelem))' recvelem];            
    [~,ind] = sortrows(dmd{i}.elemrecv); 
    %[~,ind] = sort(dmd{i}.elemrecv(:,1)); 
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
    dmd{i}.nbsend = dmd{i}.elemsend(:,1);
    dmd{i}.nbrecv = dmd{i}.elemrecv(:,1);        
    dmd{i}.elemsend = dmd{i}.elemsend(:,2);
    dmd{i}.elemrecv = dmd{i}.elemrecv(:,2);        
end

function nbelem = neighboringelements(t2t, elem)

t2te = t2t(:,elem);            
nbelem = unique(t2te(:));
if nbelem(1) == 0
    nbelem(1) = [];
end

