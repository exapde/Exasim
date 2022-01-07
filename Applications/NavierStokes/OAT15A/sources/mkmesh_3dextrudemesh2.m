function mesh = mkmesh_3dextrudemesh2(mesh2d,zz,elemtype)
% USAGE EXAMPLE: mesh = mkmesh_3dextrudemesh(mesh2d,linspace(0,0.1,10),app,elemtype);
% NOTE: It is required that mesh2d.elemtype = 1.
if mesh2d.elemtype == 0; error('mkmesh_3dextrudemesh2 can only take quad meshes as inputs'); end
if mesh2d.nodetype ~= 0; warning('mkmesh_3dextrudemesh2 may not work well for nodetype ~= 0'); end
if nargin < 3; elemtype = 1; end
porder = mesh2d.porder;
nodetype = mesh2d.nodetype;
plc1d = masternodes(porder,1,1,1);
nz = length(zz)-1;
tz = [(1:nz); (2:nz+1)]';
dz = zeros(length(plc1d),nz);
for i = 1:nz
    pz = zz(tz(i,:));
    dz(:,i) = (pz(2)-pz(1))*plc1d + pz(1);
end
if elemtype == 0 % From quad to tetra mesh
    thickness = zz(end);
    nlayers = length(zz)-1;
    [p,t,dg3d] = QtoT(mesh2d.p, mesh2d.t, thickness, nlayers, mesh2d.dgnodes, mesh2d.porder);
elseif elemtype == 1 % From quad to hexa mesh
    nxy = size(mesh2d.p,1);
    pz = repmat(zz,[nxy 1]);
    p = [repmat(mesh2d.p,[nz+1 1]) pz(:)];
    t = [];
    for i = 1:nz
        ti = [mesh2d.t+(i-1)*nxy mesh2d.t+i*nxy];
        t = [t; ti];
    end
    np2d = size(mesh2d.dgnodes,1);
    np1d = size(dz,1);
    ne2d = mesh2d.ne;
    dg3d = zeros(np2d*np1d,3,ne2d,nz);
    for i = 1:nz
        pm = repmat(dz(:,i)',[np2d 1]);
        for j = 1:ne2d
            dg3d(:,:,j,i) = [repmat(mesh2d.dgnodes(:,1:2,j),[np1d 1]) pm(:)];
        end
    end
    dg3d = reshape(dg3d,[np2d*np1d,3,ne2d*nz]);
end
disp('(p,t,dgnodes) of 3D mesh already computed...');
z1 = min(zz);
z2 = max(zz);
bndexpr = {['p(:,3)<' num2str(z1) '+1e-5'],...
           ['p(:,3)>' num2str(z2) '-1e-5'],...
           ['sqrt((p(:,1)-.5).^2+p(:,2).^2)>2'],...
           'true'};
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes = dg3d;
disp('mkmesh of 3D mesh already completed...')
end
