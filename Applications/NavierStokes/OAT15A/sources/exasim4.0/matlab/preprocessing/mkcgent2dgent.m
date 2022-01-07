function [cgnodes,cgelcon,rowent2elem,colent2elem,cgent2dgent] = mkcgent2dgent(dgnodes,tol)

dgnodes = round(dgnodes/tol)*tol;
[cgnodes,cgelcon] = mkelconcg(dgnodes);
[rowent2elem,colent2elem] = mkent2elem(cgelcon);

cgent2dgent = colent2elem;
npe = size(dgnodes,1);
nent = max(cgelcon(:));

for i = 1:nent
    nelem = rowent2elem(i+1)-rowent2elem(i);
    elem = colent2elem((rowent2elem(i)+1):rowent2elem(i+1));
    xcg = cgnodes(i,:);
    for j = 1:nelem
        xdg = dgnodes(:,:,elem(j));
        in = xiny(xcg,xdg);
        cgent2dgent(rowent2elem(i)+j) = (elem(j)-1)*npe+in;
    end
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
