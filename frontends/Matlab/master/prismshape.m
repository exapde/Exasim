function [f,fx,fy,fz] = prismshape(pts,porder)

n = size(pts,1);
if length(porder)==1
    porder = [porder porder];
end

[nf1,nf1x,nf1y]=koornwinder(pts(:,1:2),porder(1));
[g3,gz]=koornwinder(pts(:,3),porder(2)); % Legendre basis in z direction             

m = 0.5*(porder(1)+1)*(porder(1)+2);
f  = zeros(n,m*(porder(2)+1));
fx = f; fy = f; fz = f;
% perform tensor product to obtain the shape functions and their 
% derivatives on the unit cube
for ii=1:n
    f(ii,:) =  kron(g3(ii,:),nf1(ii,:));
    fx(ii,:) = kron(g3(ii,:),nf1x(ii,:));
    fy(ii,:) = kron(g3(ii,:),nf1y(ii,:));
    fz(ii,:) = kron(gz(ii,:),nf1(ii,:));
end        


