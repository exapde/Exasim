function mesh = mkmesh_ductcirc(mesh,db,H,L)
%MKMESH_DUCT Map a unit square mesh to a cos^2 Duct
%   MESH = MKMESH_DUCT(MESH,DB,DT,H)
%
%      MESH:     Mesh data structure
%                   input: mesh for the unit square created with
%                          mkmesh_square
%                   output: mesh for duct
%      DB:       Height of bottom bump
%      H:        Height of channel
%      L:       Length of channel

rb = (0.5^2+db^2)/(2*db);
x1 = 1;
x2 = 2;

a = 1;
p = L*mesh.p;
p(:,2) = loginc(p(:,2),a);
pnew(:,1) = L*p(:,1);
pnew(:,2) = H*p(:,2).*(pnew(:,1)<=x1 | pnew(:,1)>=x2) + ...
            (p(:,2).*H + ...
             (1-p(:,2)).*(-rb+db+sqrt(abs(rb^2-(pnew(:,1)-L/2).^2)))).* ...
            (pnew(:,1)>x1 & pnew(:,1)<x2);
mesh.p = pnew;

p = L*mesh.dgnodes;
p(:,2,:) = loginc(p(:,2,:),a);
pnew = zeros(size(p));
pnew(:,1,:) = L*p(:,1,:);
pnew(:,2,:) = H*p(:,2,:).*(pnew(:,1,:)<=x1 | pnew(:,1,:)>=x2) + ...
             (p(:,2,:).*(H) + ...
             (1-p(:,2,:)).*(-rb+db+sqrt(abs(rb^2-(pnew(:,1,:)-L/2).^2)))).* ...
             (pnew(:,1,:)>x1 & pnew(:,1,:)<x2);
mesh.dgnodes = pnew;
mesh.fcurved = repmat(true,size(mesh.f,1),1);
mesh.tcurved = repmat(true,size(mesh.t,1),1);


