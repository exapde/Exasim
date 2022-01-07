function [p,t] = cubemesh(m,n,o,elemtype)
%SQUAREMESH 3-D Regular TET/HEX Mesh Generator for the unit cube
%   [P,T]=SQUAREMESH(M,N,O,ELEMTYPE)
%
%      P:         Node positions (NP,2)
%      T:         Triangle indices (NT,3)
%      M:         Number of points in the x direction 
%      N:         Number of points in the y direction
%      O:         Number of points in the z direction
%      PARITY:    Flag determining the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%      ELEMTYPE:  Flag determining whether triangle or rectangle mesh
%                 Flag = 0 triangle mesh (default)
%                 Flag = 1 rectangle mesh
%   Example:
%      [p,t]=CUBEMESH(5,5,5,0);
%

if nargin<1, m=2; end
if nargin<2, n=m; end
if nargin<3, o=n; end
if nargin<4, elemtype=0; end

m = m+1;
n = n+1;
o = o+1;
    
% Generate mesh for unit cube
[x,y,z]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1),(0:o-1)/(o-1));
p=[x(:),y(:),z(:)];

if elemtype==0
    t0=[5,6,7,3; 5,6,2,3; 5,1,2,3; 6,8,7,4; 6,7,4,3; 6,2,4,3];

    map=[1,2];
    map=[map,map+m];
    map=[map,map+m*n];
    t=int32(map(t0));
    t=kron(t,ones(m-1,1,'int32'))+kron(ones(size(t),'int32'),int32(0:m-2)');
    t=kron(t,ones(n-1,1,'int32'))+kron(ones(size(t),'int32'),int32(0:n-2)'*m);
    t=kron(t,ones(o-1,1,'int32'))+kron(ones(size(t),'int32'),int32(0:o-2)'*m*n);

    %[p,t]=fixmesh(p,t);
else
    t = [1 2 m+2 m+1 m*n+1 m*n+2 m*n+m+2 m*n+m+1];
    t=kron(t,ones(o-1,1))+kron(ones(size(t)),(0:o-2)'*(m*n));        
    t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
    t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');            
end

p = p';
t = t';
