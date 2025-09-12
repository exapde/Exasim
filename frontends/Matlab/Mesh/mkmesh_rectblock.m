function mesh = mkmesh_rectblock(m,n,porder,parity,L,H,elemtype,nodetype)
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
if nargin<5, L=1;      end
if nargin<6, H=1;      end
if nargin<7, elemtype=0; end
if nargin<8, nodetype=0; end

if m < 2 || n < 2,
    error('At least m=2,n=2 needed.');
end

[p,t] = squaremesh(m,n,parity,elemtype);
p(:,1) = L*p(:,1);
p(:,2) = H*p(:,2);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all((p(:,2)>max(p0(:,2))-1e-3) & (p(:,1)>0.75*max(p0(:,1))-1e-3))',...
           'all((p(:,2)>max(p0(:,2))-1e-3) & (p(:,1)<0.75*max(p0(:,1))+1e-3) & (p(:,1)>0.25*max(p0(:,1))-1e-3))',...
           'all((p(:,2)>max(p0(:,2))-1e-3) & (p(:,1)<0.25*max(p0(:,1))+1e-3))',...
           'all(p(:,1)<1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);       

