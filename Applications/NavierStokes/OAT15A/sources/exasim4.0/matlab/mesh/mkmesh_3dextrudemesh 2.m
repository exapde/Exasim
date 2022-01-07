function mesh = mkmesh_3dextrudemesh(mesh2d,zz,hybrid)

% Ex: mesh = mkmesh_3dextrudemesh(mesh2d,linspace(0,0.1,10),'hdg');

porder   = mesh2d.porder;
elemtype = mesh2d.elemtype;
nodetype = mesh2d.nodetype;

plc1d = masternodes(porder,1,1,1);
nz = length(zz)-1;
tz = [(1:nz); (2:nz+1)]';
dz = zeros(length(plc1d),nz);
for i = 1:nz
    pz = zz(tz(i,:));
    dz(:,i) = (pz(2)-pz(1))*plc1d + pz(1);
end

nxy = size(mesh2d.p,1);
pz = repmat(zz,[nxy 1]);
p = [repmat(mesh2d.p,[nz+1 1]) pz(:)];
t = [];
for i = 1:nz
    ti = [mesh2d.t+(i-1)*nxy mesh2d.t+i*nxy];        
    t  = [t; ti]; 
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

% bndexpr = {'all(p(:,1)<min(p0(:,1))+1e-3)','all(p(:,1)>max(p0(:,1))-1e-3)', ...
%            'all(p(:,2)<min(p0(:,2))+1e-3)','all(p(:,2)>max(p0(:,2))-1e-3)', ...
%            'all(p(:,3)<min(p0(:,3))+1e-3)','all(p(:,3)>max(p0(:,3))-1e-3)'};     
%bndexpr = {'dist_NACA65_010_afterScalingAndRotation(p)<0.05','true'};
BP = 1;
alpha = (pi/180)*3.5;               % Angle of attack.
theta = (pi/180)*45 - alpha; % atan(2.5/5.0);                                                % anticlockwise rotation angle applied to airfoil
xLeftDomain = -1;                                                          % left boundary x-coord
xRightDomain = 3;                                                          % right boundary x-coord
yLeftBottom = xLeftDomain*tan(theta) - BP/2;                                                       % y-coord of bottom-left corner of blade passage mesh
yLeftTop = yLeftBottom + BP;
yRightBottom = xRightDomain*tan(theta) - BP/2;   
yRightTop = yRightBottom + BP;
z1 = min(zz);
z2 = max(zz);
bndexpr = {['p(:,3)<' num2str(z1) '+1e-3'],...
           ['p(:,3)>' num2str(z2) '-1e-3'],...
           ['p(:,1)<' num2str(xLeftDomain) '+1e-3'],...
           ['p(:,1)>' num2str(xRightDomain) '-1e-3'],...
           ['dsegment(p,[' num2str(xLeftDomain) ' ' num2str(yLeftBottom) ';' num2str(xRightDomain) ' ' num2str(yRightBottom) '])<1e-3'],...
           ['dsegment(p,[' num2str(xLeftDomain) ' ' num2str(yLeftTop) ';' num2str(xRightDomain) ' ' num2str(yRightTop) '])<1e-3'],...
           'true'};                     
mesh = mkmesh(p,t,porder,bndexpr,elemtype,nodetype);
mesh.dgnodes = dg3d;

% make periodic mesh
master = mkmaster(mesh,2*porder);
[~,mesh] = preprocess(master,mesh,hybrid);
mesh = periodic_v2(mesh,{5,'p(:,[1 3])',6,'p(:,[1 3])';1,'p(:,[1 2])',2,'p(:,[1 2])'});

% mesh.dgnodes(:,mesh2d.nd+1:2*mesh2d.nd,:) = mesh.dgnodes(:,1:mesh2d.nd,:);
% mesh.dgnodes(:,2*mesh2d.nd+1:3*mesh2d.nd,:) = 0;

end
