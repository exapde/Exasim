function mesh = periodic(mesh,bndexpr,versionFlag)

if nargin<3; versionFlag = 1; end

if strcmp(mesh.hybrid,'hdg')
    if versionFlag == 1
        mesh = periodic_v1(mesh,bndexpr);
    elseif versionFlag == 2
        mesh = periodic_v2(mesh,bndexpr);
    end
elseif strcmp(mesh.hybrid,'edg')
    if versionFlag == 1
        mesh = periodic_edgv1(mesh,bndexpr);
    elseif versionFlag == 2
        mesh = periodic_edgv2(mesh,bndexpr);
    end    
end

function mesh = periodic_v1(mesh,bndexpr)
% PERIODIC modify elcon to account for periodic boundary conditions
%
%    mesh = periodic(mesh,bndexpr)
%
%
%    bndexpr:  describes which boundaries are periodic:
%              bndexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: bndexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%          periodic, where x on bnd 1 is matched to x on bnd 2


f  = mesh.f;
t2f= mesh.t2f;

nfv= size(t2f,2);    % number of faces per element
ne = size(mesh.t,1); % number of elements
npf= size(mesh.elcon,1)/nfv; % number of points per face

[philocfc,philocvl] = localbasis(mesh.porder,mesh.nd,mesh.elemtype);
[npv,nnv] = size(philocvl);

% get dgnodes
dgnodes = zeros(npv,mesh.nd,ne);
for d=1:mesh.nd
  for n=1:nnv
    dp=philocvl(:,n)*mesh.p(mesh.t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end

elcon = mesh.elcon;
elcon = reshape(elcon,[npf nfv ne]);

for i = 1:size(bndexpr,1)    
%     ind = (mesh.f(:,end) == -bndexpr{i,1});
%     mesh.f(ind,:) = [];
%     ind = (mesh.f(:,end) == -bndexpr{i,3});
%     mesh.f(ind,:) = [];
    
    % get faces on the 1st boundary 
    [pbnd1,fbnd1,ebnd1] = getface(f,mesh.p,bndexpr{i,1},bndexpr{i,2});
    
    % get faces on the 2nd boundary 
    [pbnd2,fbnd2,ebnd2] = getface(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
    
    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1, pbnd2);    
    pbnd2 = pbnd2(ind,:);    
    fbnd2 = fbnd2(ind,:);
    ebnd2 = ebnd2(ind,:);    
    
    if max(abs(pbnd2(:)-pbnd1(:))) > 1e-10
        error('two periodic bundaries are not matched');
    end
  
    nbf = length(fbnd1);        
    % modify elcon
    for j=1:nbf
        lf1 = (t2f(ebnd1(j),:)==fbnd1(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2f(ebnd2(j),:)==fbnd2(j));  % location of face fbnd2(j) on element ebnd2(j)                
        
        tm2 = elcon(:,lf2,ebnd2(j));        
        ind = find(elcon(:)>max(tm2));
        elcon(ind) = elcon(ind) - npf;      % offset elcon to eliminate dofs on face fbnd2(j)
                
        dg2 = squeeze(dgnodes(mesh.perm(:,lf2),:,ebnd2(j)));  % dg nodes on face fbnd2(j) from element ebnd2(j)
        p   = dg2; 
        q2  = eval(bndexpr{i,4});        % evaluate periodic boundary expression   
        
        fn1 = f(fbnd1(j),1:end-2);       % corners of face fbnd1(j) 
        dg1 = philocfc*mesh.p(fn1,:);    % nodal points on face fbnd1(j)
        p   = dg1; 
        q1  = eval(bndexpr{i,2});        % evaluate periodic boundary expression   
                
        in  = xiny(q2,q1);               % match dg1 to dg2  
        tm1 = elcon(:,lf1,ebnd1(j));                
        elcon(:,lf2,ebnd2(j)) = tm1(in); % transfer dofs from face fbnd1(j) to face fbnd2(j)        
    end    
end

elcon = reshape(elcon,[npf*nfv ne]);
mesh.elcon = elcon;
mesh.nsiz = max(mesh.elcon(:));


function mesh = periodic_v2(mesh,bndexpr)
% PERIODIC modify elcon, f, t2f, nf, fcurved, bf to account for periodic boundary conditions
%
%    mesh = periodic(mesh,bndexpr)
%
%
%    bndexpr:  describes which boundaries are periodic:
%              bndexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: bndexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%
% This function only works for HDG meshes. These three lines are based on
% HDG hypotheses:
% tm2 = elcon(:,lf2,ebnd2(j));        
% ind = find(elcon(:)>max(tm2));
% elcon(ind) = elcon(ind) - npf;


tol = 1e-10; % Tolerance to match periodic boundaries.

f  = mesh.f;
t2f = mesh.t2f;
fcurved = mesh.fcurved;
bf = mesh.bf;

nfv = size(t2f,2);    % number of faces per element
ne = size(mesh.t,1); % number of elements
npf = size(mesh.elcon,1)/nfv; % number of points per face

[philocfc,philocvl] = localbasis(mesh.porder,mesh.nd,mesh.elemtype);
[npv,nnv] = size(philocvl);

% get dgnodes
dgnodes = zeros(npv,mesh.nd,ne);
for d=1:mesh.nd
  for n=1:nnv
    dp=philocvl(:,n)*mesh.p(mesh.t(:,n),d)';
    dgnodes(:,d,:)=dgnodes(:,d,:)+permute(dp,[1,3,2]);
  end
end

elcon = mesh.elcon;
elcon = reshape(elcon,[npf nfv ne]);

for i = 1:size(bndexpr,1)
%     ind = (mesh.f(:,end) == -bndexpr{i,1});
%     mesh.f(ind,:) = [];
%     ind = (mesh.f(:,end) == -bndexpr{i,3});
%     mesh.f(ind,:) = [];

    t2fUpdated = t2f;

    % get faces on the 1st boundary 
    [pbnd1,fbnd1,ebnd1] = getface(f,mesh.p,bndexpr{i,1},bndexpr{i,2});        
    
    % get faces on the 2nd boundary 
    [pbnd2,fbnd2,ebnd2] = getface(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
           
    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1, pbnd2);    
    pbnd2 = pbnd2(ind,:);    
    fbnd2 = fbnd2(ind,:);
    ebnd2 = ebnd2(ind,:);    
    
    if max(abs(pbnd2(:)-pbnd1(:))) > tol
        error('two periodic bundaries are not matched');
    end
  
    nbf = length(fbnd1);        

    for j=1:nbf
        disp(['j= ',num2str(j),' / ',num2str(nbf)]);
        
        lf1 = (t2f(ebnd1(j),:)==fbnd1(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2f(ebnd2(j),:)==fbnd2(j));  % location of face fbnd2(j) on element ebnd2(j)                
        
        % Update elcon:
        tm2 = elcon(:,lf2,ebnd2(j));        
        ind = find(elcon(:)>max(tm2));
        elcon(ind) = elcon(ind) - npf;      % offset elcon to eliminate dofs on face fbnd2(j)

        p   = squeeze(dgnodes(mesh.perm(:,lf1),:,ebnd1(j)));  % dg nodes on face fbnd1(j) from element ebnd1(j)        
        q1  = eval(bndexpr{i,2});        % evaluate periodic boundary expression           
%         fn1 = f(fbnd1(j),1:end-2);       % corners of face fbnd1(j) 
%         p   = philocfc*mesh.p(fn1,:);    % nodal points on face fbnd1(j)        
%         q1  = eval(bndexpr{i,2});        % evaluate periodic boundary expression   
        
        p   = squeeze(dgnodes(mesh.perm(:,lf2),:,ebnd2(j)));  % dg nodes on face fbnd2(j) from element ebnd2(j)        
        q2  = eval(bndexpr{i,4});        % evaluate periodic boundary expression   
        
        in  = xiny(q2,q1);               % match q2 to q1  
        tm1 = elcon(:,lf1,ebnd1(j));     % dofs of face fbnd1(j)           
        elcon(:,lf2,ebnd2(j)) = tm1(in); % transfer dofs from face fbnd1(j) to face fbnd2(j)
        
        % Update f:
        f(fbnd1(j),end) = f(fbnd2(j),end-1);
        % f(fbnd1(j),end) = ebnd2(j);
        % f(fbnd2(j),end) = ebnd1(j);
        
        % Update t2fUpdated:
        t2fUpdated(ebnd2(j),lf2) = t2fUpdated(ebnd1(j),lf1);
        % t2fUpdated(ebnd2(j),lf2) = fbnd1(j);
        
        bf(lf1,ebnd1(j)) = ebnd2(j);
        bf(lf2,ebnd2(j)) = ebnd1(j);
        
    end
    
    % Decrease entries in t2fUpdated to account for the faces removed:            
    facesToBeRemoved = sort(fbnd2,'descend');    % Faces to be removed, ordered from largest to smallest index.
    facesToBeRemoved = reshape(facesToBeRemoved,1,[]);
    
    for ff = facesToBeRemoved
        if abs(max(max(t2fUpdated == ff))) ~= 0
            error('A face that needed to be removed was not properly removed from t2f.');
        end
        t2fUpdated = t2fUpdated - (t2fUpdated > ff);
    end
    
    % Remove rows in f and fcurved:
    f = f(setdiff(1:size(f,1),fbnd2),:);
    fcurved = fcurved(setdiff(1:length(fcurved),fbnd2),:);
    
    t2f = t2fUpdated;
    
end

mesh.elcon = reshape(elcon,[npf*nfv ne]);
mesh.f = f;
mesh.t2f = t2f;
mesh.f2f = mkf2f(f, t2f);
mesh.nf = size(mesh.f,1);
mesh.fcurved = fcurved;
mesh.nsiz = max(mesh.elcon(:));
mesh.bf = bf;


function [pbnd,fbnd,ebnd] = getface(f,p0,bn,expr)

% get faces on the boundary bn
fbnd = find(f(:,end)==-bn); 
ft   = f(fbnd,:);
% number of faces on the boundary bn
n    = length(fbnd); 

p = p0(ft(1,1:end-2),:); 
q = eval(expr); 
[k, m] = size(q); % m = nd, k = number of vertices of a face

% get elements which contain faces on the boundary bn
ebnd = ft(:,end-1);

% get coordinate values of faces on the boundary bn
pbnd = zeros(k,m,n);
for i = 1:n % for each face i on the boundary bn
    p = p0(ft(i,1:end-2),:); % vertices of the face i
    q = eval(expr); % evaluate periodic expression of this face   
    pbnd(:,:,i) = sortrows(snap(q));    
end
pbnd = reshape(permute(pbnd,[3,1,2]),n,[]);

function x=snap(x)

tol=sqrt(eps);
x=tol*round(x/tol);


function in = xiny(x,y)
% Determine if each row of x is a member of y
% If row j of x is a member of y and x(j,:) = y(k,:) then in(j) = k
% Else in(j) = 0

[m,dim] = size(x);
in = zeros(m,1);
if dim==1
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
elseif dim==2
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
elseif dim==3
    for j=1:m
        d2 = (y(:,1)-x(j,1)).^2 + (y(:,2)-x(j,2)).^2 + (y(:,3)-x(j,3)).^2;
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end
else
    n = size(y,1);    
    for j=1:m
        d2 = sum((y - repmat(x(j,:),[n 1])).^2,2);
        [md,id] = min(d2);
        if md<1e-12, in(j)=id; end
    end    
end

function [philocfc,philocvl] = localbasis(porder,dim,elemtype) 

[plocvl,tlocal,plocfc] = masternodes(porder,dim,elemtype,0);

if dim==2 && elemtype==0      % tri
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = 1 - xi - eta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
elseif dim==2 && elemtype==1  % quad
    xi  = plocfc(:,1);
    philocfc(:,1) = 1 - xi;
    philocfc(:,2) = xi;
    xi  = plocvl(:,1);
    eta = plocvl(:,2);    
    philocvl(:,1) = (1-xi).*(1-eta);
    philocvl(:,2) = xi.*(1-eta);
    philocvl(:,3) = xi.*eta;
    philocvl(:,4) = (1-xi).*eta;
elseif dim==3 && elemtype==0  % tet
    xi  = plocfc(:,1);
    eta = plocfc(:,2);    
    philocfc(:,1) = 1 - xi - eta;
    philocfc(:,2) = xi;
    philocfc(:,3) = eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = 1 - xi - eta - zeta;
    philocvl(:,2) = xi;
    philocvl(:,3) = eta;
    philocvl(:,4) = zeta;
elseif dim==3 && elemtype==1   % hex
    xi  = plocfc(:,1);
    eta = plocfc(:,2);
    philocfc(:,1) = (1-xi).*(1-eta);
    philocfc(:,2) = xi.*(1-eta);
    philocfc(:,3) = xi.*eta;
    philocfc(:,4) = (1-xi).*eta;
    xi   = plocvl(:,1);
    eta  = plocvl(:,2);
    zeta = plocvl(:,3);
    philocvl(:,1) = (1-xi).*(1-eta).*(1-zeta);
    philocvl(:,2) = xi.*(1-eta).*(1-zeta);
    philocvl(:,3) = xi.*eta.*(1-zeta);
    philocvl(:,4) = (1-xi).*eta.*(1-zeta);    
    philocvl(:,5) = (1-xi).*(1-eta).*(zeta);
    philocvl(:,6) = xi.*(1-eta).*(zeta);
    philocvl(:,7) = xi.*eta.*(zeta);
    philocvl(:,8) = (1-xi).*eta.*(zeta);        
end

function mesh = periodic_edgv1(mesh,bndexpr)
% PERIODIC modify elcon, f, t2f, nf, fcurved, bf to account for periodic boundary conditions
%
%    mesh = periodic(mesh,bndexpr)
%
%
%    bndexpr:  describes which boundaries are periodic:
%              bndexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: bndexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%

tol = 1e-10; % Tolerance to match periodic boundaries.

f  = mesh.f;
t2f = mesh.t2f;

nfv = size(t2f,2);    % number of faces per element
ne = size(mesh.t,1); % number of elements
npf = size(mesh.elcon,1)/nfv; % number of points per face

% get elcon and edgnodes
mesh.f(:,end+1) = 1;
[elcon,~,edgnodes] = elconnectivities(mesh);
elcon = reshape(elcon,[npf nfv ne]);
mesh.f(:,end) = [];

for i = 1:size(bndexpr,1)

    % get faces on the 1st boundary 
    [pbnd1,fbnd1,ebnd1] = getface(f,mesh.p,bndexpr{i,1},bndexpr{i,2});             
    %[pbnd1,fbnd1,ebnd1] = getedgn(f,mesh.p,bndexpr{i,1},bndexpr{i,2});             
    
    % get faces on the 2nd boundary 
    [pbnd2,fbnd2,ebnd2] = getface(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
    %[pbnd2,fbnd2,ebnd2,fand2] = getedgn(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
    
    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1, pbnd2);    
    pbnd2 = pbnd2(ind,:);    
    fbnd2 = fbnd2(ind,:);    
    ebnd2 = ebnd2(ind,:);        
    
    if max(abs(pbnd2(:)-pbnd1(:))) > tol
        error('two periodic bundaries are not matched');
    end
             
    % update elcon, t2f, bf   
    for j=1:length(ebnd2)
        lf1 = (t2f(ebnd1(j),:)==fbnd1(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2f(ebnd2(j),:)==fbnd2(j));  % location of face fbnd2(j) on element ebnd2(j)                        
        i1 = elcon(:,lf1,ebnd1(j));
        i2 = elcon(:,lf2,ebnd2(j));
        p = edgnodes(i1,:);
        q1 = eval(bndexpr{i,2}); % evaluate periodic boundary expression           
        p = edgnodes(i2,:);
        q2 = eval(bndexpr{i,4}); % evaluate periodic boundary expression   
        in = xiny(q2,q1);   % match q2 to q1          
        elcon(:,lf2,ebnd2(j)) = elcon(in,lf1,ebnd1(j));
    end    
    
    % edgnodes to be removed, ordered from largest to smallest index.    
    edgNodesToBeRemoved = [];
    for j=1:length(ebnd2)
        edgNodesToBeRemoved = [edgNodesToBeRemoved; elcon(:,fand2(j),ebnd2(j))];
    end    
    edgNodesToBeRemoved = unique(edgNodesToBeRemoved);
    edgNodesToBeRemoved = sort(edgNodesToBeRemoved,'descend');       
    for j=1:length(edgNodesToBeRemoved)
        %ind = find(elcon(:)>edgNodesToBeRemoved(j));
        %elcon(ind) = elcon(ind) - 1;
        elcon = elcon - (elcon > edgNodesToBeRemoved(j));
    end            
end
mesh.elcon = reshape(elcon,[npf*nfv ne]);
mesh.nsiz = max(mesh.elcon(:));


function mesh = periodic_edgv2(mesh,bndexpr)
% PERIODIC modify elcon, f, t2f, nf, fcurved, bf to account for periodic boundary conditions
%
%    mesh = periodic(mesh,bndexpr)
%
%
%    bndexpr:  describes which boundaries are periodic:
%              bndexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: bndexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%

tol = 1e-10; % Tolerance to match periodic boundaries.

f  = mesh.f;
t2f = mesh.t2f;
fcurved = mesh.fcurved;
bf = mesh.bf;

nfv = size(t2f,2);    % number of faces per element
ne = size(mesh.t,1); % number of elements
npf = size(mesh.elcon,1)/nfv; % number of points per face

% get elcon and edgnodes
mesh.f(:,end+1) = 1;
[elcon,~,edgnodes] = elconnectivities(mesh);
elcon = reshape(elcon,[npf nfv ne]);
mesh.f(:,end) = [];

for i = 1:size(bndexpr,1)

    % get faces on the 1st boundary 
    [pbnd1,fbnd1,ebnd1] = getface(f,mesh.p,bndexpr{i,1},bndexpr{i,2});             
    %[pbnd1,fbnd1,ebnd1] = getedgn(f,mesh.p,bndexpr{i,1},bndexpr{i,2});             
    
    % get faces on the 2nd boundary 
    [pbnd2,fbnd2,ebnd2] = getface(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
    %[pbnd2,fbnd2,ebnd2,fand2] = getedgn(f,mesh.p,bndexpr{i,3},bndexpr{i,4});        
    
    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1, pbnd2);    
    pbnd2 = pbnd2(ind,:);    
    fbnd2 = fbnd2(ind,:);    
    ebnd2 = ebnd2(ind,:);        
    
    if max(abs(pbnd2(:)-pbnd1(:))) > tol
        error('two periodic bundaries are not matched');
    end
          
    t2fUpdated = t2f;
    elconUpdated = elcon;
    
    % update elcon, t2f, bf, f, fcurved
    edgNodesToBeRemoved = [];
    for j=1:length(ebnd2)
        disp(['j= ',num2str(j),' / ',num2str(length(ebnd2))]);
        
        lf1 = (t2f(ebnd1(j),:)==fbnd1(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2f(ebnd2(j),:)==fbnd2(j));  % location of face fbnd2(j) on element ebnd2(j)                        
        i1 = elcon(:,lf1,ebnd1(j));
        i2 = elcon(:,lf2,ebnd2(j));
        newEdgNodesToBeRemoved = elcon(:,lf2,ebnd2(j));
        edgNodesToBeRemoved = [edgNodesToBeRemoved; newEdgNodesToBeRemoved(:)];
        p = edgnodes(i1,:);
        q1 = eval(bndexpr{i,2}); % evaluate periodic boundary expression           
        p = edgnodes(i2,:);
        q2 = eval(bndexpr{i,4}); % evaluate periodic boundary expression   
        in = xiny(q2,q1);   % match q2 to q1
        
        % Update elcon
        for k=1:length(newEdgNodesToBeRemoved)
            indices = find(elcon(:) == newEdgNodesToBeRemoved(k));
            elconUpdated(indices) = elcon(in(k),lf1,ebnd1(j));
        end
        % elcon(:,lf2,ebnd2(j)) = elcon(in,lf1,ebnd1(j));
        
        % Update f:
        f(fbnd1(j),end) = f(fbnd2(j),end-1);
        
        % Update t2fUpdated:
        t2fUpdated(ebnd2(j),lf2) = t2fUpdated(ebnd1(j),lf1);
        
        bf(lf1,ebnd1(j)) = ebnd2(j);
        bf(lf2,ebnd2(j)) = ebnd1(j);
    end    
    
    % edgnodes to be removed, ordered from largest to smallest index.    
    edgNodesToBeRemoved = unique(edgNodesToBeRemoved);
    edgNodesToBeRemoved = sort(edgNodesToBeRemoved,'descend');       
    for j=1:length(edgNodesToBeRemoved)
        edgnodes = [edgnodes(1:edgNodesToBeRemoved(j)-1,:);edgnodes(edgNodesToBeRemoved(j)+1:end,:)];
        elconUpdated = elconUpdated - (elconUpdated > edgNodesToBeRemoved(j));
    end    
    
    % Faces to be removed, ordered from largest to smallest index.    
    facesToBeRemoved = sort(fbnd2,'descend');       
    for j=1:length(facesToBeRemoved)
        %ind = find(t2f(:)>facesToBeRemoved(j));
        %t2f(ind) = t2f(ind) - 1;
        t2fUpdated = t2fUpdated - (t2fUpdated > facesToBeRemoved(j));
    end     
    
    f = f(setdiff(1:size(f,1),fbnd2),:);
    fcurved = fcurved(setdiff(1:length(fcurved),fbnd2),:);
    
    t2f = t2fUpdated;
    elcon = elconUpdated;
end

mesh.elcon = reshape(elcon,[npf*nfv ne]);
mesh.f = f;
mesh.t2f = t2f;
mesh.f2f = mkf2f(f, t2f);
mesh.nf = size(mesh.f,1);
mesh.fcurved = fcurved;
mesh.nsiz = max(mesh.elcon(:));
mesh.bf = bf;


% function [edgn,fbnd,ebnd] = getedgn(mesh,bn)
% 
% % get all edgnodes 
% mesh.f(:,end+1)=1;
% [elcon,~,edg] = elconnectivities(mesh);
% mesh.f = [];
% 
% % get faces on the boundary bn
% fbnd = find(mesh.f(:,end)==-bn); 
% ft   = mesh.f(fbnd,:);
% % number of faces on the boundary bn
% n    = length(fbnd); 
% 
% % get elements which contain faces on the boundary bn
% ebnd = ft(:,end-1);
% 
% % find edgnodes on periodic boundaries
% edgn = [];
% for i = 1:length(ebnd)
%     el = reshape(elcon(:,ebnd(i)),size(mesh.perm));
%     t2f = mesh.t2f(i,:);
%     for j = 1:length(t2f)
%         if ismember(t2f(j),ft)
%             edgn = [edgn; edg(el(:,j),:)];
%         end
%     end
% end
% edgn = sortrows(snap(edgn));  
