function [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(dgnodes,tol,opts)

if nargin<3
    opts = 0;
end

dgnodes = round(dgnodes/tol)*tol;

% CG element-to-node connectivities
[cgnodes,cgelcon] = mkelconcg(dgnodes);

% CG node-to-element connectivities
[rowent2elem,colent2elem] = mkent2elem(cgelcon);

cgent2dgent = colent2elem;
npe = size(dgnodes,1);
nent = max(cgelcon(:));

for i = 1:nent % for CG node i
    nelem = rowent2elem(i+1)-rowent2elem(i); % number of elements connected to CG node i    
    elem = colent2elem((rowent2elem(i)+1):rowent2elem(i+1)); % elements connected to CG node i
    xcg = cgnodes(i,:); % coordinates of CG node i
    for j = 1:nelem % for element j connected to CG node i
        xdg = dgnodes(:,:,elem(j)); % coordinates of DG nodes on element j
        in = xiny(xcg,xdg); % match CG node i to one of the DG nodes on element j
        cgent2dgent(rowent2elem(i)+j) = (elem(j)-1)*npe+in; % index of the matched DG node
    end
end

if (opts==1)
% Modify [rowent2elem,colent2elem] so that all CG nodes are connected to
% the same number of elements
colent2elem = mkcgent2elem(dgnodes, cgnodes, cgelcon, colent2elem, rowent2elem);
rowent2elem = size(colent2elem,1)*ones(nent,1);
rowent2elem = [0; cumsum(rowent2elem)];
colent2elem = colent2elem(:);
end

% 
% ne = size(cgelcon,2);
% for m = 1:ne
%     tm = cgnodes(cgelcon(:,m),:)-dgnodes(:,:,m);
%     if max(abs(tm(:)))>tol
%         error('mkelconcg is wrong');
%     end
% end
% 
% nd = size(cgnodes,2);
% dgnodes = reshape(permute(dgnodes,[1 3 2]),[npe*ne nd]);
% for i = 1:nent
% %     nelem = rowent2elem(i+1)-rowent2elem(i);
% %     elem = colent2elem((rowent2elem(i)+1):rowent2elem(i+1));
%     dgn = cgent2dgent((rowent2elem(i)+1):rowent2elem(i+1));
%     xcg = cgnodes(i,:);
%     xdg = dgnodes(dgn,:);
% %     xdg = zeros(nelem,nd);
% %     for j = 1:nelem
% %         xdg(j,:) = dgnodes(dgn(j),:,elem(j));
% %     end
%     tm = xcg - xdg;
%     if max(abs(tm(:)))>tol
%         error('mkcgent2dgent is wrong');
%     end
% end
% 
% 
% 
% 
% 
