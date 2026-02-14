function mesh = mkmesh_sphericalfrustum(slant, porder, n, L, alpha, dz)
% slant: coordinates of an inclined line of the frustum

% a quarter of the unit circle mesh
mesh2d = mkmesh_quartercircle(porder,n);

if nargin < 4 || isempty(L)
  points = sphericalfrustum(mesh2d.p, slant);
else
  points = sphericalfrustum(mesh2d.p, slant, L, alpha);
end

pn = squeeze(points(:,:,1));
tn = [];
np = size(pn,2);
ns = size(slant,2);
t0 = mesh2d.t;
for i = 2:ns
    pn = cat(2, pn, squeeze(points(:,:,i)));
    t1 = t0 + np;
    tloc = cat(1, t0, t1);
    tn = cat(2, tn, tloc);
    t0 = t1;
end

mesh = mkmesh(pn',tn',porder,{'true'},1,1);
mesh.p = mesh.p';
mesh.t = mesh.t';

nz = ns-1;
if nargin < 6
  plc1d = masternodes(porder,1,1);
  tz = [(1:nz); (2:nz+1)]';
  dz = zeros(3,length(plc1d),nz);
  for i = 1:nz
      pz = slant(:,tz(i,:));
      dz(:,:,i) = (pz(:,2)-pz(:,1))*plc1d' + pz(:,1);
  end
end
dz = reshape(dz, 3, (porder+1)*nz);

[npe2d, nd2d, ne2d] = size(mesh2d.dgnodes);
p = reshape(permute(mesh2d.dgnodes, [2 1 3]), nd2d, npe2d*ne2d);

if nargin < 4 || isempty(L)
  points = sphericalfrustum(p, dz);
else
  points = sphericalfrustum(p, dz, L, alpha);
end

points = reshape(points, [3 npe2d ne2d porder+1 nz]);
points = permute(points, [2 4 1 3 5]);
mesh.dgnodes = reshape(points, [npe2d*(porder+1) 3 ne2d*nz]);

figure(1); clf; meshplot(mesh,1);
axis on;

end

function points = sphericalfrustum(p, slant, L, alpha)

% the radii of the cross sections of the frustum
R = abs(slant(2,:));

% the radius of the left base
R1 = R(1);    

% the radius of the right base
R2 = R(end);  

% the heigh of the frustum
H = abs(slant(1,1) - slant(1,end)); 

if nargin < 3
  % the height of the inner cone formed by the left base
  L = H/(R2/R1-1); 
end

if nargin < 4  
  alpha = 0; 
end

ns = size(slant,2);
np = size(p,2);
points = zeros(3, np, ns);
for i = 1:ns
  r = R(i);        %  radius of the conical frustum cross section
  x = slant(1,i);  %  x-coordinate of the conical frustum cross section  
  y = r*p(2,:);    %  y-coordinate of points on the conical frustum cross section   
  z = r*p(1,:);    %  z-coordinate of points on the conical frustum cross section   
  a = sqrt(y.^2 + z.^2);   % radius of the original circle 
  d = abs(x - slant(1,1)); % distance from the top base of radius r
  h = L+alpha*d;                 % the height of the cone whose base is radius r
  l = sqrt(r^2 + h^2);     % slant height of the cone whose base is radius r
  m = sqrt(a.^2 + h^2);    % slant height of the inner cone whose base is radius a
  b = (a*l)./m;            % radius of the projected circle is computed by triangle similarity 
  s = (h*l)./m;            % slant height of the projected cone whose base is radius b
  points(2,:,i) = b.*p(2,:); % y-coordinate of points on the projected circle of radius b
  points(3,:,i) = b.*p(1,:); % z-coordinate of points on the projected circle of radius b
  points(1,:,i) = x - (s - h); % x-coordinate of points on the projected circle of radius b
end

end