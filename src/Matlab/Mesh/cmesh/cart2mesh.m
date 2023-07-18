function mesh = cart2mesh(porder,X,Y,Z,bndexpr,elemtype,nodetype)
% Create DG mesh from a Cartesian grid

if nargin<7; nodetype=0; end

if isempty(Z) 
  [p,t,p1]=cart2dg(elemtype,porder,X,Y);
else
  [p,t,p1]=cart2dg(elemtype,porder,X,Y,Z);
end

mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes = p1;

