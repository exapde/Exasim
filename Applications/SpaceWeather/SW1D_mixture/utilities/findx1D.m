function [indx,xlocal] = findx1D(x,dgnodes)

nelems = size(dgnodes,3);

dgnodesminmax = reshape(dgnodes([1,end],1,:),[2,nelems])';

[indx,~] = find((x'>dgnodesminmax(:,1)).*(x'<dgnodesminmax(:,2)));
xlocal = (x-dgnodesminmax(indx,1))./(dgnodesminmax(indx,2)-dgnodesminmax(indx,1));