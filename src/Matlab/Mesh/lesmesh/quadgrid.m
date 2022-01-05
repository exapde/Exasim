function [p,t] = quadgrid(x,y)

[X,Y] = ndgrid(x,y);
p = [X(:) Y(:)];

m = length(x);
n = length(y);
t = [1 2 m+2 m+1];
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');

% m = length(xv)-1;
% n = length(yv)-1;
% t = [1, 2, m+3, m+2];
% t = kron(t,ones(m,1)) + kron(ones(size(t)),(0:m-1)');
% t = kron(t,ones(n,1)) + kron(ones(size(t)),(0:n-1)'*(m+1));

