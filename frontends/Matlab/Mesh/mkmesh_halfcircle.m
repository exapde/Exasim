function mesh = mkmesh_halfcircle(mesh,a,b,c,t1,t2)
%MKMESH_DUCT Map a unit square mesh to a cos^2 Duct
%   MESH = MKMESH_DUCT(MESH,DB,DT,H)
%
%      MESH:     Mesh data structure
%                   input: mesh for the unit square created with
%                          mkmesh_square
%                   output: mesh for duct
%      DB:       Height of bottom bump
%      DT:       Height of top bump
%      H:        Height of channel
%      L:        Length of channel
%
p = mesh.p;
t = t2 + p(:,2)*(t1-t2);
d = -b*cos(t) + c*(1+cos(t));
r = d + p(:,1).*(a-d);

pnew(:,1) = r.*cos(t);
pnew(:,2) = r.*sin(t);
mesh.p = pnew;

if isfield(mesh,'dgnodes') && ~isempty(mesh.dgnodes)
   clear pnew;
   p = mesh.dgnodes;
   t = t2 + p(:,2,:)*(t1-t2);
   d = -b*cos(t) + c*(1+cos(t));
   r = d + p(:,1,:).*(a-d);   
   pnew = zeros(size(p));
   pnew(:,1,:) = r.*cos(t);
   pnew(:,2,:) = r.*sin(t);
   mesh.dgnodes = pnew;
   mesh.fcurved = true(size(mesh.f,1),1);
   mesh.tcurved = true(size(mesh.t,1),1);
end

