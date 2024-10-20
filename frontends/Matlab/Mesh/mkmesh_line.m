function mesh = mkmesh_line(m,porder,a,b,nodetype)
%MKMESH_LINE Creates 1D mesh data structure for in the interval [a,b].
%   MESH=MQMESH_LINE(M,PORDER,a,b,nodetype)
%
%      MESH:      Mesh structure
%      M:         Number of points 
%      PORDER:    Polynomial Order of Approximation (default=1)
%      NODETYPE:  Flag determining node distribution 
%                 Flag = 0 uniform distribution (default)
%                 Flag = 1 nonuniform distribution

%   See also: MKMESH
%

if nargin<1, m=2; end
if nargin<2, porder=1; end
if nargin<3, a=0;      end
if nargin<4, b=1;      end
if nargin<5, nodetype=0;  end

[p,t]=masternodes(m-1,1);
p = a + (b-a)*p;

bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)'};     
elemtype = 0;

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);

