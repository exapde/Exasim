function [p,t,dgnodes] = sphereCube(order, R0, R1, n, m)

[pc,tc] = cubemesh(n,n,m,1);

pc(3,:) = loginc(pc(3,:),3);

pc(1,:) = -pi/4 + pi*pc(1,:)/2;
pc(2,:) = -pi/4 + pi*pc(2,:)/2;
pc(3,:) = R0 + (R1-R0)*pc(3,:);

dgnod = createdgnodes(pc,tc,0,[],[],order);
[npe,nd,ne] = size(dgnod);
dgnod = reshape(permute(dgnod,[2,1,3]),nd,npe*ne);

offe = size(tc,2);
offp = size(pc,2);
offd = npe*ne;

ph = zeros(3,6*offp);
th = zeros(8,6*offe);
dg = zeros(3,6*npe*ne);

p = transform(pc);
g = transform(dgnod);

ph(:,       1:  offp) = p;
dg(:,       1:  offd) = g;
th(:,       1:  offe) = tc;

% figure
% meshplot3D(ph,tc,1)
% hold on


R = rotmat([0,0,1]);
p2 = R*p;
g2 = R*g;
ph(:,  offp+1:2*offp) = p2;
dg(:,  offd+1:2*offd) = g2;
th(:,  offe+1:2*offe) = tc +   offp;
% meshplot3D(ph,tc + offp,1)
% hold on
p3 = R*p2;
g3 = R*g2;
ph(:,2*offp+1:3*offp) = p3;
dg(:,2*offd+1:3*offd) = g3;
th(:,2*offe+1:3*offe) = tc + 2*offp;
p4 = R*p3;
g4 = R*g3;
ph(:,3*offp+1:4*offp) = p4;
dg(:,3*offd+1:4*offd) = g4;
th(:,3*offe+1:4*offe) = tc + 3*offp;

R = rotmat([0,1,0]);
p5 = R'*p;
g5 = R'*g;
ph(:,4*offp+1:5*offp) = p5;
dg(:,4*offd+1:5*offd) = g5;
th(:,4*offe+1:5*offe) = tc + 4*offp;
p6 = R*p;
g6 = R*g;
ph(:,5*offp+1:6*offp) = p6;
dg(:,5*offd+1:6*offd) = g6;
th(:,5*offe+1:6*offe) = tc + 5*offp;

% meshplot3D(ph,tc + 5*offp,1)

[ph,th] = fixmesh(ph',th');

p = ph';
t = th';

dgnodes = permute(reshape(dg,nd,npe,6*ne),[2,1,3]);

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

function pp = transform(p)
tanth = tan(p(1,:));
tanph = tan(p(2,:));
x = p(3,:)./sqrt(1+tanth.^2+tanph.^2);
pp = [x; x.*tanth; x.*tanph];
end

function y=loginc(x,alpha)
a=min(x(:));b=max(x(:));
y = a + (b-a)*(exp(alpha.*(x-a)/(b-a))-1)./(exp(alpha)-1);
end
