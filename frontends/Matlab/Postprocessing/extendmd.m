function dmd = extendmd(dmd)

% number of subdomains
nproc = length(dmd);

% faces on the interface of subdomain i   
intf = cell(nproc,1);         
for i = 1:nproc
    ne1 = dmd{i}.elempartpts(1) + dmd{i}.elempartpts(2); % number of elements within subdomain i   
    ind1 = find(dmd{i}.f2t(1,:) > ne1); 
    ind2 = find(dmd{i}.f2t(3,:) > ne1);
    intf{i} = unique([ind1 ind2]); % faces on the interface of subdomain i   
end

% neighboring subdomains of subdomain i for each external element received 
nbsdelemrecv = cell(nproc,1); 
for i = 1:nproc
  nbsdelemrecv{i} = 0*dmd{i}.elemrecv;
  n = 0;
  for j = 1:length(dmd{i}.nbsd)    
    ind = (n+1):(n+dmd{i}.elemrecvpts(j));
    nbsdelemrecv{i}(ind) = dmd{i}.nbsd(j);    
    n = n + dmd{i}.elemrecvpts(j); 
  end
end

for i = 1:nproc
  dmd{i}.interfaces = zeros(length(intf{i}), 16);
  for j = 1:length(intf{i}) % loop over each face on the interface of subdomain i
    fj = dmd{i}.f2t(:,intf{i}(j));     
    e1 = fj(1); % 1st local element
    e2 = fj(3); % 2nd local element
    g1 = dmd{i}.elempart(e1); % 1st global element on subdomain i
    g2 = dmd{i}.elempart(e2); % 2nd global element on subdomain i

    ee = max(e1, e2);         % external element
    ii = dmd{i}.elemrecv==ee; % find ee in dmd{i}.elemrecv
    n = nbsdelemrecv{i}(ii);  % n is a neighboring subdomain of i

    for k = 1:length(intf{n}) % loop over each face on the interface of subdomain n
      fk = dmd{n}.f2t(:,intf{n}(k));     
      h1 = dmd{n}.elempart(fk(1)); % 1st global element on subdomain n
      h2 = dmd{n}.elempart(fk(3)); % 2nd global element on subdomain n       
      if ((g1==h1) && (g2==h2)) || ((g1==h2) && (g2==h1))
        % collect interface information when the global elements match
        dmd{i}.interfaces(j,:) = [i intf{i}(j) fj' g1 g2 n intf{n}(k) fk' h1 h2];     
        break;
      end        
    end            
  end
end
