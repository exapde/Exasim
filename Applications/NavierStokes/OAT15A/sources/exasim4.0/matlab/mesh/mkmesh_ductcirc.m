function mesh = mkmesh_ductcirc(mesh,db,H)
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

rb = (0.5^2+db^2)/(2*db);

p = mesh.p;
pnew(:,1) = 3*p(:,1);
pnew(:,2) = H*p(:,2).*(pnew(:,1)<=1 | pnew(:,1)>=2) + ...
            (p(:,2).*H + ...
             (1-p(:,2)).*(-rb+db+sqrt(rb^2-(pnew(:,1)-1.5).^2))).* ...
            (pnew(:,1)>1 & pnew(:,1)<2);
mesh.p = pnew;

if isfield(mesh,'dgnodes') & ~isempty(mesh.dgnodes)
   clear pnew;
   p = mesh.dgnodes;
   pnew = zeros(size(p));
   pnew(:,1,:) = 3*p(:,1,:);
   pnew(:,2,:) = H*p(:,2,:).*(pnew(:,1,:)<=1 | pnew(:,1,:)>=2) + ...
                 (p(:,2,:).*(H) + ...
                 (1-p(:,2,:)).*(-rb+db+sqrt(rb^2-(pnew(:,1,:)-1.5).^2))).* ...
                 (pnew(:,1,:)>1 & pnew(:,1,:)<2);
   mesh.dgnodes = pnew;
   mesh.fcurved = repmat(true,size(mesh.f,1),1);
   mesh.tcurved = repmat(true,size(mesh.t,1),1);
end
