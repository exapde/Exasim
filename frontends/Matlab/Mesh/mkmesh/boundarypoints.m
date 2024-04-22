function pb = boundarypoints(p,f,ib)

dim = size(p,2);
ii = f(:,end)==-ib;
idx = f(ii,1:dim);
idx = unique(idx(:));
pb = p(idx,:);  
