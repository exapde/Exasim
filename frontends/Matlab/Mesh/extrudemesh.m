function mesh = extrudemesh(mesh2d,zz)

porder = mesh2d.porder;
plc1d = masternodes(porder,1,1);

nz = length(zz)-1;
tz = [(1:nz); (2:nz+1)]';
dz = zeros(length(plc1d),nz);
for i = 1:nz
    pz = zz(tz(i,:));
    dz(:,i) = (pz(2)-pz(1))*plc1d + pz(1);
end

nxy = size(mesh2d.p',1);
pz = repmat(zz,[nxy 1]);
mesh.p = [repmat(mesh2d.p',[nz+1 1]) pz(:)]';

t2d = mesh2d.t';
[ne2d, nv2d] = size(t2d);

mesh.t = zeros(ne2d*nz, nv2d*2);
for i = 1:nz
    mesh.t(ne2d*(i-1)+1:ne2d*i,:) = [t2d+(i-1)*nxy t2d+i*nxy];    
end
mesh.t = mesh.t';

np2d = size(mesh2d.dgnodes,1);
np1d = size(dz,1);
mesh.dgnodes = zeros(np2d*np1d,3,ne2d,nz);
for i = 1:nz
    pm = repmat(dz(:,i)',[np2d 1]);
    for j = 1:ne2d
        mesh.dgnodes(:,:,j,i) = [repmat(mesh2d.dgnodes(:,1:2,j),[np1d 1]) pm(:)];
    end
end
mesh.dgnodes = reshape(mesh.dgnodes,[np2d*np1d,3,ne2d*nz]);


