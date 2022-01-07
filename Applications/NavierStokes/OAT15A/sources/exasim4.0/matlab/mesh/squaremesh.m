function [p,t]=squaremesh(m,n,parity,elemtype)
%SQUAREMESH 2-D Regular Triangle/Rectangle Mesh Generator for the unit square
%   [P,T]=SQUAREMESH(M,N,PARITY,ELEMTYPE)
%
%      P:         Node positions (NP,2)
%      T:         Triangle indices (NT,3)
%      PARITY:    Flag determining the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%      ELEMTYPE:  Flag determining whether triangle or rectangle mesh
%                 Flag = 0 triangle mesh (default)
%                 Flag = 1 rectangle mesh
%   Example:
%      [p,t]=SQUAREMESH(5,10,1);
%

if nargin<1, m=10; end
if nargin<2, n=m; end
if nargin<3, parity=0; end
if nargin<4, elemtype=0; end

% Generate mesh for unit square

[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1));
p=[x(:),y(:)];

if elemtype==0 % triangular mesh
    if parity==0
      t=[1,2,m+2; 1,m+2,m+1];
    else
      t=[1,2,m+1; 2,m+2,m+1];
    end

    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);

    % Reorder triangles in Cartesian order

    ix=[];
    for i=1:n-1
      ix1=i+(n-1)*(0:m-2);
      ix2=ix1+(n-1)*(m-1);

      if parity==0
        ix12=[ix2;ix1];
      else
        ix12=[ix1;ix2];
      end

      ix=[ix,reshape(ix12,1,[])];
    end

    t=t(ix,:);
else % rectangular mesh
    t = [1 2 m+2 m+1];
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
end
