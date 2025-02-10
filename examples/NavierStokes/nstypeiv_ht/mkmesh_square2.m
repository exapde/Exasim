function mesh = mkmesh_square2(m,n,porder,parity,L,H,elemtype,nodetype)
%MKMESH_SQUARE Creates 2D tri/quad mesh data structure for unit square.
%   MESH=MKMESH_SQUARE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the horizaontal direction 
%      N:         Number of points in the vertical direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      PARITY:    Flag determining the the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%      ELEMTYPE:  Flag determining element type
%                 Flag = 0 tri/tet elements (default)
%                 Flag = 1 quad/hex elements
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution

%   See also: SQUAREMESH, MKMESH
%

if nargin<2, m=2; n=2; end
if nargin<3, porder=1; end
if nargin<4, parity=0; end
if nargin<5, L=1;      end
if nargin<6, H=1;      end
if nargin<7, elemtype=0; end
if nargin<8, nodetype=1; end

%elemtype = 1;
if m < 2 || n < 2,
    error('At least m=2,n=2 needed.');
end

[p,t] = squaremesh2(m,n,parity,elemtype);
p(:,1) = L*p(:,1);
p(:,2) = H*p(:,2);

bndexpr = {'all(p(:,2)<1e-6)','all(p(:,1)>max(p0(:,1))-1e-6)', ...
           'all(p(:,2)>max(p0(:,2))-1e-6)','all(p(:,1)<1e-6)'}; 

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.p = p';
mesh.t = t';

% mesh = mkmesh_halfcircle(mesh, 1, 3, 4.7, pi/2, 3*pi/2);
% mesh.porder = porder;
% mesh.boundaryexpr = {@(p) sqrt(p(1,:).^2+p(2,:).^2)<1+1e-6, @(p) p(1,:)>-1e-7, @(p) abs(p(1,:))<20};
% mesh.bndexpr = {'all(sqrt(p(:,1).^2+p(:,2).^2)<1+1e-6)'    'all(p(:,1)>-1e-7)'    'all(abs(p(:,1))<20)'}
% mesh.periodicexpr = {};

function [p,t]=squaremesh2(m,n,parity,elemtype)
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
x = logdec(x,4);
x1 = unique(x(:));
[~,j] = min(abs(x1-0.95));
ii = x<=x1(j); x(ii) = loginc(x(ii),2);
ii = x>=x1(j); x(ii) = logdec(x(ii),2.0);

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

