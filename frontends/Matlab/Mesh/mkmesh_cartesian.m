function mesh = mkmesh_cartesian(x,y,porder,parity,elemtype,nodetype)
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

[p,t] = rectmesh(x,y,parity,elemtype);
p = p';
t = t';

bndexpr = {'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
           'all(p(:,2)>max(p0(:,2))-1e-3)','all(p(:,1)<min(p0(:,1))+1e-3)'};     

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.p = p';
mesh.t = t';


