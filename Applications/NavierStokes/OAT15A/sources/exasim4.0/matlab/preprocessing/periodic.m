function [tp,fp,t2fp,t2tp,faceconp] = periodic(p,t,f,t2f,t2t,facecon,periodicexpr)
% PERIODIC modify facecon to account for periodic boundary conditions
%
%    mesh = periodic(mesh,periodicexpr)
%
%
%    periodicexpr:  describes which boundaries are periodic:
%              periodicexpr = {bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: periodicexpr = {1,'p(:,1),3,'p(:,1)'} will make bnds 1,3
%          periodic, where x on bnd 1 is matched to x on bnd 3

tp = t;
fp = f;
t2fp = t2f;
t2tp = t2t;
faceconp = facecon;

nperiodic = size(periodicexpr,1);
pbnd1 = cell(nperiodic,1);
fbnd1 = cell(nperiodic,1);
ebnd1 = cell(nperiodic,1);
vbnd1 = cell(nperiodic,1);
xbnd1 = cell(nperiodic,1);
pbnd2 = cell(nperiodic,1);
fbnd2 = cell(nperiodic,1);
ebnd2 = cell(nperiodic,1);
vbnd2 = cell(nperiodic,1);
xbnd2 = cell(nperiodic,1);
for i = 1:nperiodic
    % get nodes on the 1st boundary 
    [pbnd1{i},vbnd1{i},fbnd1{i},ebnd1{i},xbnd1{i}] = getface(f,p,periodicexpr{i,1},periodicexpr{i,2});        
    
    % get nodes on the 2nd boundary 
    [pbnd2{i},vbnd2{i},fbnd2{i},ebnd2{i},xbnd2{i}] = getface(f,p,periodicexpr{i,3},periodicexpr{i,4});                       
    
    % sort pbnd2 to match it to pbnd1
    ind = xiny(pbnd1{i}, pbnd2{i});    
    pbnd2{i} = pbnd2{i}(ind,:);    
    vbnd2{i} = vbnd2{i}(ind,:);      
    
    % sort xbnd2 to match it to zbnd1
    ind = xiny(xbnd1{i}, xbnd2{i});    
    xbnd2{i} = xbnd2{i}(ind,:);    
    fbnd2{i} = fbnd2{i}(ind,:);
    ebnd2{i} = ebnd2{i}(ind,:);        
        
    if max(abs(pbnd2{i}(:)-pbnd1{i}(:))) > 1e-10
        error('two periodic bundaries are not matched');
    end
    if max(abs(xbnd2{i}(:)-xbnd1{i}(:))) > 1e-10
        error('two periodic bundaries are not matched');
    end    
end

% eleminate DOFs on the 2nd boundary and keep DOFs on the 1st boundary
for i = 1:size(periodicexpr,1)       
    
    % face-to-element (indices)
    fp(fbnd1{i},end) = fp(fbnd2{i},end-1);   
    
    % face-to-element (DG nodes)
    faceconp(:,2,fbnd1{i}) = faceconp(:,1,fbnd2{i});
    
    % set DOFs on the 2nd boundary to the DOFs on the 1st boundary
    for j = 1:length(vbnd2{i})                
        tp(tp(:)==vbnd2{i}(j)) = vbnd1{i}(j);        
    end    
    
    nbf = length(fbnd1{i});        
    for j=1:nbf        
        lf1 = (t2f(ebnd1{i}(j),:)==fbnd1{i}(j));  % location of face fbnd1(j) on element ebnd1(j)
        lf2 = (t2f(ebnd2{i}(j),:)==fbnd2{i}(j));  % location of face fbnd2(j) on element ebnd2(j)                        
        t2fp(ebnd2{i}(j),lf2) = t2fp(ebnd1{i}(j),lf1);    
        t2tp(ebnd2{i}(j),lf2) = ebnd1{i}(j);
        t2tp(ebnd1{i}(j),lf1) = ebnd2{i}(j);    
    end    
end
fp(cell2mat(fbnd2),:) = [];
faceconp(:,:,cell2mat(fbnd2)) = [];

mapping = zeros(max(tp(:)),1);
uniqueNodes = unique(tp(:));
mapping(uniqueNodes) = 1:length(uniqueNodes);
tp = mapping(tp);

mapping = zeros(max(t2fp(:)),1);
uniqueNodes = unique(t2fp(:));
mapping(uniqueNodes) = 1:length(uniqueNodes);
t2fp = mapping(t2fp);

function [pbnd,vbnd,fbnd,ebnd,xbnd] = getface(f,p0,bn,expr)

% get faces on the boundary bn
fbnd = find(f(:,end)==-bn); 
% get elements on the boundary bn
ebnd = f(fbnd,end-1);
% get vertices on the boundary bn
vbnd = f(fbnd,1:end-2);
vbnd = unique(vbnd(:));
% get nodes on the boundary bn
p = p0(vbnd,:); 
pbnd = eval(expr); 
pbnd = snap(pbnd);

% number of faces on the boundary bn
n = length(fbnd); 
ft   = f(fbnd,:);
p = p0(ft(1,1:end-2),:); 
q = eval(expr); 
[k, m] = size(q); % k = number of vertices of a face
% get coordinate values of faces on the boundary bn
xbnd = zeros(k,m,n);
for i = 1:n % for each face i on the boundary bn
    p = p0(ft(i,1:end-2),:); % vertices of the face i
    q = eval(expr); % evaluate periodic expression of this face   
    xbnd(:,:,i) = sortrows(snap(q));    
end
xbnd = reshape(permute(xbnd,[3,1,2]),n,[]);


