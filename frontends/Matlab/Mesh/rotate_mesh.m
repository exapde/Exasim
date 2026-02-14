function [p, t, xdg] = rotate_mesh(p0, t0, ns, xdg0)

dth = 0.5*pi/ns;
%p0 = mesh.p;
p0 = cat(1, p0, p0(1,:)*0);

%t0 = mesh.t;
np = size(p0,2);

p = p0;
t = [];
th = 0;
for i = 1:ns
    th = th + dth;
    R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
    p = cat(2, p, R*p0);
    t1 = t0 + np;
    tn = cat(1, t0, t1);
    t = cat(2, t, tn);
    t0 = t1;
end

if nargin > 3
  [npe2d, nd2d, ne2d] = size(xdg0);
  p0 = reshape(permute(xdg0, [2 1 3]), [nd2d npe2d*ne2d]);
  p0 = cat(1, p0, p0(1,:)*0);

  porder = sqrt(npe2d) - 1;
  plc1d = masternodes(porder,1,1);
  ta = linspace(0, 0.5*pi, ns+1);
  ts = [(1:ns); (2:ns+1)]';
  theta = zeros(length(plc1d),ns);
  for i = 1:ns
      pz = ta(ts(i,:));
      theta(:,i) = (pz(2)-pz(1))*plc1d + pz(1);
  end  
  theta = theta(:);

  xdg = p0;
  for i = 2:length(theta)
      th = theta(i);
      R = [1 0 0; 0 cos(th) -sin(th); 0 sin(th) cos(th)];
      xdg = cat(2, xdg, R*p0);
  end  
  xdg = reshape(xdg, [3 npe2d ne2d (porder+1) ns]);
  xdg = permute(xdg, [2 4 1 3 5]);
  xdg = reshape(xdg, [npe2d*(porder+1) 3 ne2d*ns]);
else
  xdg = [];
end

end