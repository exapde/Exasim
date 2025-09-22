function [dmd1, dmd2, isd1, isd2] = interfacepartition(mesh1, dmd1, mesh2, dmd2)

[isd1, xdg1, xm1] = interfacesubdomain(mesh1, dmd1);
[isd2, xdg2, xm2] = interfacesubdomain(mesh2, dmd2);

[dmd1, dmd2, match, nbsd] = interfacematching(xm1, xm2, dmd1, dmd2, isd1, isd2);

n1 = length(xm1); % number of interface subdomains in Domain 1
n2 = length(xm2); % number of interface subdomains in Domain 2
l = 1;
for i = 1:n1
  for j = 1:n2
    if nbsd(i,j)==1
      k1 = match{i,j}(:,1);
      k2 = match{i,j}(:,2);      
      x1 = xdg1{i}(:,:,k1);
      x2 = xdg2{j}(:,:,k2);
      for m = 1:length(k1)
        if (l==1)
          perm1 = xiny(x1(:,:,m), x2(:,:,m));
          perm2 = xiny(x2(:,:,m), x1(:,:,m));        
        else
          tm1 = xiny(x1(:,:,m), x2(:,:,m));
          tm2 = xiny(x2(:,:,m), x1(:,:,m));      
          if (max(abs(perm1-tm1)) > 0) || (max(abs(perm2-tm2)) > 0)
            error("something wrong");
          end
        end
        l = l + 1;
      end
    end
  end
end

for i = 1:n1
  d = isd1(i);
  dmd1{d}.faceperm = perm1;
  dmd1{d}.nbintf = dmd1{d}.nbintf + length(dmd1); 
end

for i = 1:n2
  d = isd2(i);
  dmd2{d}.faceperm = perm2;
end

end

function [isd, xdg, xm] = interfacesubdomain(mesh, dmd)

coupledinterface = mesh.coupledinterface;
porder = mesh.porder;
perm = mesh.perm;
npf = size(perm,1);
nd = size(mesh.p,1);

% find subdomains containing the coupled interface
isd = zeros(length(dmd),1); 
for i = 1:length(dmd)
  if (dmd{i}.intepartpts(2) > 0)
    isd(i) = 1;
  end  
end
isd = find(isd==1); % list of interface subdomains

% for each interface subdomain, find all DG nodes on the interface
xdg = cell(length(isd),1);
xm = cell(length(isd),1);
for i = 1:length(isd)
  n = isd(i); % subdomain n
  m = (dmd{n}.intepartpts(1)+1):(dmd{n}.intepartpts(1)+dmd{n}.intepartpts(2));
  e = dmd{n}.elempart(m); % a list of all the interface elements of subdomain n 
  
  % DG nodes for interface elements 
  dgnodes = createdgnodes(mesh.p, mesh.t(:,e), mesh.f(:,e), [], [], porder);        
  
  % find local faces on the interface
  [intf, ~] = find(mesh.f(:,e) == coupledinterface);
    
  ne = length(e);    % number of faces on the interface
  xdg{i} = zeros(npf, nd, ne);  % DG nodes for faces on the interface
  for j = 1:ne
    xdg{i}(:,:,j) = dgnodes(perm(:,intf(j)),:,j);
  end    
  
  % the center points for faces on the interface
  xm{i} = reshape(mean(xdg{1},1), [nd ne])';    
end

end

function [dmd1, dmd2, match, nbsd] = interfacematching(xm1, xm2, dmd1, dmd2, isd1, isd2)

n1 = length(xm1); % number of interface subdomains in Domain 1
n2 = length(xm2); % number of interface subdomains in Domain 2

match = cell(n1, n2);
nbsd = zeros(n1, n2);
for i = 1:n1    % for each interface subdomain i in Domain 1
  x1 = xm1{i};  % interface nodes on subdomain i in Domain 1
  for j = 1:n2  % for each interface subdomain j in Domain 2
    x2 = xm2{j};% interface nodes on subdomain j in Domain 2
    in = xiny(x1, x2); % match interface nodes
    k1 = find(in > 0); % faces belong to subdomain i 
    k2 = in(k1);       % faces belong to subdomain j
    if isempty(k1)     
      match{i,j} = []; % subdomain i and subdomain j are not connected
    else
      match{i,j} = [k1(:) k2(:)]; % collect interface faces on subdomain i and subdomain j
      nbsd(i,j) = 1; % subdomain i and subdomain j are connected
    end
  end  
end

nbintf1 = cell(n1 , 1);      % store neigboring subdomains
facesend1 = cell(n1 , 1);    % store faces sent to neigboring subdomains
facesendpts1 = cell(n1 , 1); % store numbers of faces sent to neigboring subdomains
facerecv1 = cell(n1 , 1);    % store faces received from neigboring subdomains
facerecvpts1 = cell(n1 , 1); % store numbers of faces received from neigboring subdomains
for i = 1:n1
  facesend1{i} = [];
  facerecv1{i} = [];
  nbintf1{i} = find(nbsd(i,:) == 1);  
  nnb = length(nbintf1{i});
  facesendpts1{i} = zeros(nnb, 1);
  facerecvpts1{i} = zeros(nnb, 1);
end

nbintf2 = cell(n2, 1);       % store neigboring subdomains
facesend2 = cell(n2 , 1);    % store faces sent to neigboring subdomains
facesendpts2 = cell(n2 , 1); % store numbers of faces sent to neigboring subdomains
facerecv2 = cell(n2 , 1);    % store faces received from neigboring subdomains
facerecvpts2 = cell(n2 , 1); % store numbers of faces received from neigboring subdomains
for j = 1:n2
  facesend2{j} = [];
  facerecv2{j} = [];
  nbintf2{j} = find(nbsd(:,j) == 1);
  nnb = length(nbintf2{j});
  facesendpts2{j} = zeros(nnb, 1);
  facerecvpts2{j} = zeros(nnb, 1);
end

for i = 1:n1
  nnb = length(nbintf1{i});
  for k = 1:nnb
    j = nbintf1{i}(k);
    m = length(match{i,j}(:,1));
    facesend1{i} = [facesend1{i}; [j*ones(m,1) match{i,j}(:,1)]];    
    facerecv2{j} = [facerecv2{j}; [i*ones(m,1) match{i,j}(:,2)]];
    facesendpts1{i}(k) = m;
    facerecvpts2{j}(nbintf2{j}==i) = m;
  end
end

for j = 1:n2
  nnb = length(nbintf2{j});
  for k = 1:nnb
    i = nbintf2{j}(k);    
    m = length(match{i,j}(:,1));    
    facesend2{j} = [facesend2{j}; [i*ones(m,1) match{i,j}(:,2)]];    
    facerecv1{i} = [facerecv1{i}; [j*ones(m,1) match{i,j}(:,1)]];
    facesendpts2{j}(k) = m;
    facerecvpts1{j}(nbintf1{i}==j) = m;
  end
end

for i = 1:n1
  d = isd1(i);
  dmd1{d}.nbintf = isd2(nbintf1{i});
  dmd1{d}.facesend = [isd2(facesend1{i}(:,1)) facesend1{i}(:,2)];
  dmd1{d}.facesendpts = facesendpts1{i};
  dmd1{d}.facerecv = [isd2(facerecv1{i}(:,1)) facerecv1{i}(:,2)];
  dmd1{d}.facerecvpts = facerecvpts1{i};
end

for i = 1:n2
  d = isd2(i);
  dmd2{d}.nbintf = isd1(nbintf2{i});
  dmd2{d}.facesend = [isd1(facesend2{i}(:,1)) facesend2{i}(:,2)];
  dmd2{d}.facesendpts = facesendpts2{i};
  dmd2{d}.facerecv = [isd1(facerecv2{i}(:,1)) facerecv2{i}(:,2)];
  dmd2{d}.facerecvpts = facerecvpts2{i};
end

end
