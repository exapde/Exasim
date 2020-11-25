function [p,t] = sphere_cube_uniform(R0, R1, n, m)

[pc,tc] = cubemesh(n,n,m,1);

pc(1,:) = -pi/4 + pi*pc(1,:)/2;
pc(2,:) = -pi/4 + pi*pc(2,:)/2;
pc(3,:) = R0 + (R1-R0)*pc(3,:);

offe = size(tc,2);
offp = size(pc,2);

ph = zeros(3,6*offp);
th = zeros(8,6*offe);

p = trasnform(pc);

ph(:,       1:  offp) = p;
th(:,       1:  offe) = tc;

R = rotmat([0,0,1]);
p2 = R*p;
ph(:,  offp+1:2*offp) = p2;
th(:,  offe+1:2*offe) = tc +   offp;
p3 = R*p2;
ph(:,2*offp+1:3*offp) = p3;
th(:,2*offe+1:3*offe) = tc + 2*offp;
p4 = R*p3;
ph(:,3*offp+1:4*offp) = p4;
th(:,3*offe+1:4*offe) = tc + 3*offp;

R = rotmat([0,1,0]);
p5 = R'*p;
ph(:,4*offp+1:5*offp) = p5;
th(:,4*offe+1:5*offe) = tc + 4*offp;
p6 = R*p;
ph(:,5*offp+1:6*offp) = p6;
th(:,5*offe+1:6*offe) = tc + 5*offp;

[ph,th] = fixmesh(ph',th');

p = ph';
t = th';

end

function R = rotmat(u);
W = [    0  -u(3)   u(2);
      u(3)      0  -u(1);
     -u(2)   u(1)     0];

R = eye(3,3) + W + W^2;
end

function [p,t]=fixmesh(p,t)
%FIXMESH  Remove duplicated/unused nodes 
% Remove duplicated nodes
snap = 100*max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx] = unique(round(p/snap)*snap,'rows');
p = p(ix,:);
t = jx(t);
if size(t,2) == 1, t = t'; end  % This lines ensures the function works for one element

% Remove nodes that are not contained in t:
[pix,ix,jx] = unique(t);
t = reshape(jx,size(t));
p = p(pix,:);
end

function pp = trasnform(p)
tanth = tan(p(1,:));
tanph = tan(p(2,:));
x = p(3,:)./sqrt(1+tanth.^2+tanph.^2);
pp = [x; x.*tanth; x.*tanph];
end
