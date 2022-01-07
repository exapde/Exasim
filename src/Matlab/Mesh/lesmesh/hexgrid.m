function [p,t] = hexgrid(x,y,z)

[X,Y,Z]=ndgrid(x,y,z);
p=[X(:),Y(:),Z(:)];

m = length(x);
n = length(y);
o = length(z);

t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
t=kron(t,ones(o-1,1))+kron(ones(size(t)),(0:o-2)'*(m*n));        
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');            

