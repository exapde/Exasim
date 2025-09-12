function mesh = mkmesh_rect(m,n,porder,parity,xrect,elemtype,nodetype)
%MKMESH_SQUARE Creates 2D mesh data structure for unit square.
%   MESH=MKMESH_SQUARE(M,N,PORDER,PARITY)
%
%      MESH:      Mesh structure
%      M:         Number of points in the horizaontal direction 
%      N:         Number of points in the vertical direction
%      PORDER:    Polynomial Order of Approximation (default=1)
%      PARITY:    Flag determining the the triangular pattern
%                 Flag = 0 (diagonals SW - NE) (default)
%                 Flag = 1 (diagonals NW - SE)
%
%   See also: SQUAREMESH, MKT2F, SETBNDNBRS, UNIFORMLOCALPNTS, CREATENODES
%

if nargin<2, m=2; n=2; end
if nargin<3, porder=1; end
if nargin<4, parity=0; end
if nargin<5, xrect=[0 1 0 1];  end

if m < 2 || n < 2,
    error('At least m=2,n=2 needed.');
end

[p,t] = squaremesh(m,n,parity,elemtype);
p(:,1) = xrect(1) + (xrect(2)-xrect(1))*p(:,1);
p(:,2) = xrect(3) + (xrect(4)-xrect(3))*p(:,2);

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
