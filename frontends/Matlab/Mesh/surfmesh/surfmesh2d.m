function mesh = surfmesh2d(xl, xu, ns, np, porder, ss, sp, alpha, dgs, dgp)
%SURFMESH creates a two-dimensional DG mesh defined by the two curves
%        __________________________________  
%       |                                  | xu 
%       |                                  |
%    np |                                  |
%       |                                  |
%       |__________________________________| xl
%                       ns
%
%   xl   : points defines the bottom curve
%   xu   : points defines the top curve
%   ns   : number of subdivisons along the curves  
%   np   : number of subdivisions along the normal direction 
% porder : polynomial degree
%   ss   : scaling factors along the curves  
%   sp   : scaling factors along the normal direction 
%  alpha : parameter to control the grid distribution along the normal direction

% Below are the steps to compute p, t, and dgnodes 
% step 1  : get 2 end points from xl -> xl1 and xl2
% step 2  : get 2 end points from xu -> xu1 and xu2
% step 3  : map [xl1, xl2] to the unit interval
% step 3  : map [xu1, xu2] to the unit interval
% step 4  : map xl to the bottom edge of the unit square
% step 4  : map xu to the uppwer edge of the unit square

% step 3 : map these 4 end points to a unit square

if nargin<5
  porder = 1;
  ss = [0 0];
  sp = [0 0];
  alpha = 0;
  dgs = [];
  dgp = [];
elseif nargin<6
  ss = [0 0];
  sp = [0 0];
  alpha = 0;
  dgs = [];
  dgp = [];
elseif nargin<8
  alpha = 0;
  dgs = [];
  dgp = [];
elseif nargin<9  
  dgs = [];
  dgp = [];  
end

% make sure ss and sp positive 
for i = 1:2
  ss(i) = max(ss(i), 1e-8);
  sp(i) = max(sp(i), 1e-8);
end

if numel(ns)>1
  xs = ns;
  ns = size(xs,1)-1;
else
  xs = [];
end

if numel(np)>1
  xp = np;
  np = size(xp,1)-1;
else
  xp = [];
end

% figure(1); clf; plot(xl(:,1), xl(:,2), 'o', 'LineWidth', 2);
% hold on;
% plot(xu(:,1), xu(:,2), 's', 'LineWidth', 2);

nl = size(xl, 1);
nu = size(xu, 1);

% 4 end points
xl1 = xl(1,:);    % bottom left point
xl2 = xl(end,:);  % bottom right point
xu1 = xu(1,:);    % top left point
xu2 = xu(end,:);  % top right point

% map these 4 end points to a unit square
% 2-D bilinear map: xl1 -> (0,0), xl2 -> (1,0), xu2 -> (1, 1), xu1 -> (0, 1)    
% x = c(1,1) + c(2,1) * xi + c(3,1) * eta + c(4,1) * xi * eta
% y = c(1,2) + c(2,2) * xi + c(3,2) * eta + c(4,2) * xi * eta
A = [1 0 0 0; 1 1 0 0; 1 1 1 1; 1 0 1 0];
b = [xl1; xl2; xu2; xu1];
c = A\b;

% find points in the reference domain that map to xl
zl = 0*xl;
for i = 1:nl
  [zl(i,1), zl(i,2)] = inverse_bilinear_mapping(xl(i,1), xl(i,2), c, (i-1)/(nl-1), 0);
end

% use polyfit to find the polynomial that fits zl
Pl = polyfit(zl(:,1) , zl(:,2) , min(nl-1, 12));

% find points in the reference domain that map to xu
zu = 0*xu;
for i = 1:nu
  [zu(i,1), zu(i,2)] = inverse_bilinear_mapping(xu(i,1), xu(i,2), c, (i-1)/(nu-1), 1);
end

if ~isempty(xs)
  zs = 0*xs;
  n = size(zs,1);
  for i = 1:n
    [zs(i,1), zs(i,2)] = inverse_bilinear_mapping(xs(i,1), xs(i,2), c, (i-1)/(n-1), 0);
  end  
end

if ~isempty(dgs)  
  n = size(dgs,1);
  for i = 1:n
    [dgs(i,1), dgs(i,2)] = inverse_bilinear_mapping(dgs(i,1), dgs(i,2), c, (i-1)/(n-1), 0);
  end  
end

if ~isempty(xp)
  zp = 0*xp;
  n = size(zp,1);
  for i = 1:n
    [zp(i,1), zp(i,2)] = inverse_bilinear_mapping(xp(i,1), xp(i,2), c, 0, (i-1)/(n-1));
  end  
end

if ~isempty(dgp)  
  n = size(dgp,1);
  for i = 1:n
    [dgp(i,1), dgp(i,2)] = inverse_bilinear_mapping(dgp(i,1), dgp(i,2), c, 0, (i-1)/(n-1));
  end  
end

% use polyfit to find the polynomial that fits zu
Pu = polyfit(zu(:,1) , zu(:,2) , min(nu-1, 12));

if isempty(xs) && isempty(xp)
  % construct a uniform grid on the unit square
  mesh = mkmesh_square(ns, np, porder, 1, 1, 1, 1, 1);
elseif isempty(xp)
  mesh = mkmesh_cartesian(zs(:,1), linspace(0,1,np+1)', porder, 1, 1, 1);
elseif isempty(xs)
  mesh = mkmesh_cartesian(linspace(0,1,ns+1)', zp(:,2), porder, 1, 1, 1);  
else
  mesh = mkmesh_cartesian(zs(:,1), zp(:,2), porder, 1, 1, 1);
end

if ~isempty(dgs)
  % map dgs to mesh.dgnodes(:,1,:)  
  dgs = reshape(dgs, porder+1, ns, 2);
  dg = mesh.dgnodes(:,1,:);
  dg = reshape(dg, [porder+1, porder+1, ns, np]);
  for n = 1:ns
    for i = 1:(porder+1)
      dg(i,:,n,:) = dgs(i,n,1);
    end
  end
  mesh.dgnodes(:,1,:) = reshape(dg, [(porder+1)*(porder+1) 1 ns*np]);
end

if ~isempty(dgp)
  % map dgp to mesh.dgnodes(:,2,:)  
  dgp = reshape(dgp, porder+1, np, 2);
  dg = mesh.dgnodes(:,2,:);
  dg = reshape(dg, [porder+1, porder+1, ns, np]);
  for n = 1:np
    for j = 1:(porder+1)
      dg(:,j,:,n) = dgp(j,n,2);
    end
  end
  mesh.dgnodes(:,2,:) = reshape(dg, [(porder+1)*(porder+1) 1 ns*np]);  
end

% reshape the grid points
cg0 = mesh.p';
[npe, dim, ne] = size(mesh.dgnodes);
dg0 = zeros(npe*ne, dim);
dg0(:,1) = reshape(mesh.dgnodes(:,1,:), [npe*ne 1]);
dg0(:,2) = reshape(mesh.dgnodes(:,2,:), [npe*ne 1]);

bottom = abs(cg0(:,2) - 0) < 1e-6;
top = abs(cg0(:,2) - 1) < 1e-6;
left = abs(cg0(:,1) - 0) < 1e-6;
right = abs(cg0(:,1) - 1) < 1e-6;

dbottom = abs(dg0(:,2) - 0) < 1e-6;
dtop = abs(dg0(:,2) - 1) < 1e-6;
dleft = abs(dg0(:,1) - 0) < 1e-6;
dright = abs(dg0(:,1) - 1) < 1e-6;


% modify the grid distribution along the eta direction
if alpha > 0
  % cg0(:,2) = (1- alpha*cg0(:,1)) .* cg0(:,2) + alpha*cg0(:,1) .* ((1 - cos(pi * cg0(:,2))) / 2);
  % dg0(:,2) = (1- alpha*dg0(:,1)) .* dg0(:,2) + alpha*dg0(:,1) .* ((1 - cos(pi * dg0(:,2))) / 2);

  cg0(:,2) = (1- alpha*cg0(:,1)) .* cg0(:,2) + alpha*cg0(:,1) .* logdec(cg0(:,2),7);
  dg0(:,2) = (1- alpha*dg0(:,1)) .* dg0(:,2) + alpha*dg0(:,1) .* logdec(dg0(:,2),7);
end

% apply log scaling distribution to xi direction
cg0(:,1) = logdec(loginc(cg0(:,1), ss(1)), ss(2));
dg0(:,1) = logdec(loginc(dg0(:,1), ss(1)), ss(2));

% apply log scaling distribution to eta direction
cg0(:,2) = logdec(loginc(cg0(:,2), sp(1)), sp(2));
dg0(:,2) = logdec(loginc(dg0(:,2), sp(1)), sp(2));

% map eta to y using y = yl + (yu - yl) * eta
yl = polyval(Pl, cg0(:,1)); % compute yl using the polynomial of the bottom curve
yu = polyval(Pu, cg0(:,1)); % compute yu using the polynomial of the uppwer curve
ind = abs(cg0(:,1))<1e-6 | abs(cg0(:,1)-1)<1e-6; 
yl(ind) = 0;
yu(ind) = 1;
cg0(:,2) = yl + (yu - yl).*cg0(:,2);

yl = polyval(Pl, dg0(:,1,:));
yu = polyval(Pu, dg0(:,1,:));
ind = abs(dg0(:,1))<1e-6 | abs(dg0(:,1)-1)<1e-6; 
yl(ind) = 0;
yu(ind) = 1;
dg0(:,2) = yl + (yu - yl).*dg0(:,2);

% figure(1); clf; plot(zl(:,1), zl(:,2), 'o', 'LineWidth', 2);
% hold on;
% plot(zu(:,1), zu(:,2), 's', 'LineWidth', 2);
% plot(cg0(:,1), cg0(:,2), 'x');
% axis equal; axis tight;

% Apply forward bilinear mapping 
[p(:,1), p(:,2)] = forward_bilinear_mapping(cg0(:,1), cg0(:,2), c);
[dgnodes(:,1), dgnodes(:,2)] = forward_bilinear_mapping(dg0(:,1), dg0(:,2), c);

mesh.cgbottom = p(bottom,:);
mesh.cgtop = p(top,:);
mesh.cgleft = p(left,:);
mesh.cgright = p(right,:);

mesh.dgbottom = dgnodes(dbottom,:);
mesh.dgtop = dgnodes(dtop,:);
mesh.dgleft = dgnodes(dleft,:);
mesh.dgright = dgnodes(dright,:);

mesh.cbottom = bottom;
mesh.ctop = top;
mesh.cleft = left;
mesh.cright = right;
mesh.dbottom = dbottom;
mesh.dtop = dtop;
mesh.dleft = dleft;
mesh.dright = dright;

% figure(2); clf; plot(xl(:,1), xl(:,2), 'o', 'LineWidth', 2);
% hold on;
% plot(xu(:,1), xu(:,2), 's', 'LineWidth', 2);
% plot(p(:,1), p(:,2), 'x');
% axis equal; axis tight;

p = p';
t = mesh.t;
dgnodes = reshape(dgnodes, [npe, ne, dim]);
dgnodes = permute(dgnodes, [1 3 2]);

mesh.xpe = mesh.plocal;
mesh.p = p;
mesh.t = t;
mesh.dgnodes = dgnodes;

% figure(3); clf; meshplot(mesh);
% hold on;
% plot(xl(:,1), xl(:,2), 'o', 'LineWidth', 2);
% plot(xu(:,1), xu(:,2), 's', 'LineWidth', 2);
% axis equal; axis tight;

end

% t = linspace(pi, pi/4, 100)';
% xl = cos(t);
% yl = sin(t);
% xu = 2*cos(t);
% yu = 2*sin(t);
% [p, t, dgnodes, mesh] = surfmesh2d([xl yl], [xu yu], 30, 30, 3);



% 
% % xy = [1 xi(:) eta(:) xi(:).*eta(:)]*c;
% 
% % C = x(:,1)'/[1,1,1,1; xi1; xi2; xi1.*xi2];
% % D = x(:,2)'/[1,1,1,1; xi1; xi2; xi1.*xi2];
%  
% px=p(:,1); py=p(:,2);
% p=[0*px+1, px, py, px.*py]*[C;D]';
% 
% 
% 
% zl = 0*xl;
% for i = 1:nl
%   [x1, x2, position] = find_perpendicular_point(xl1(1), xl1(2), xl2(1), xl2(2), xl(i,1), xl(i,2));
%   zl(i,1) = sqrt((x1 - xl1(1))^2 + (x2 - xl1(2))^2);
%   zl(i,2) = position * sqrt((x1 - xl(i,1))^2 + (x2 - xl(i,2))^2);    
% end
% %dl = sqrt((xl2(1) - xl1(1))^2 + (xl2(2) - xl1(2))^2);
% 
% zu = 0*xu;
% for i = 1:nu
%   [x1, x2, position] = find_perpendicular_point(xu1(1), xu1(2), xu2(1), xu2(2), xu(i,1), xu(i,2));
%   zu(i,1) = sqrt((x1 - xu1(1))^2 + (x2 - xu1(2))^2);
%   zu(i,2) = position * sqrt((x1 - xu(i,1))^2 + (x2 - xu(i,2))^2);    
% end
% %du = sqrt((xu2(1) - xu1(1))^2 + (xu2(2) - xu1(2))^2);
% 
% 
% 
% Ns = ns*porder + 1;
% Np = np*porder + 1;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
