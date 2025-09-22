function plocal = tetnodes(p)
% TRINODES compute optimized nodal points on the master tetrahedron
%      PLOCAL = TETNODES(PORDER)
%
%      PORDER:    Order of the complete polynomial                  
%
%      PLOCAL:    Node positions (NPL,2)
%                 NPL = (PORDER+1)*(PORDER+2)*(PORDER+3)/6

% choose optimized blending parameter
alphastore = [0;0;0;0.1002; 1.1332;1.5608;1.3413;1.2577;1.1603;...
              1.10153;0.6080;0.4523;0.8856;0.8717;0.9655];
if(p<=15); alpha = alphastore(p) ; else  alpha = 1. ; end 

% total number of nodes and tolerance
%N = (p+1)*(p+2)*(p+3)/6; 
tol = 1e-10;
[r,s,t] = EquiNodes3D(p); % create equidistributed nodes
L1 = (1+t)/2; L2 = (1+s)/2; L3 = -(1+r+s+t)/2; L4 = (1+r)/2;

% set vertices of tetrahedron
v1 = [-1, -1/sqrt(3), -1/sqrt(6)]; v2 = [ 1, -1/sqrt(3),-1/sqrt(6)]; 
v3 = [ 0,  2/sqrt(3), -1/sqrt(6)]; v4 = [ 0,  0,         3/sqrt(6)];

% orthogonal axis tangents on faces 1-4
t1(1,:) = v2-v1;          t1(2,:) = v2-v1;
t1(3,:) = v3-v2;          t1(4,:) = v3-v1;   
t2(1,:) = v3-0.5*(v1+v2); t2(2,:) = v4-0.5*(v1+v2);
t2(3,:) = v4-0.5*(v2+v3); t2(4,:) = v4-0.5*(v1+v3);  

for n=1:4 % normalize tangents
   t1(n,:) = t1(n,:)/norm(t1(n,:)); t2(n,:) = t2(n,:)/norm(t2(n,:)); 
end  


% Warp and blend for each face (accumulated in shiftXYZ)
XYZ = L3*v1+L4*v2+L2*v3+L1*v4; % form undeformed coordinates
shift = zeros(size(XYZ));
for face=1:4 
  if(face==1); La = L1; Lb = L2; Lc = L3; Ld = L4; end;
  if(face==2); La = L2; Lb = L1; Lc = L3; Ld = L4; end;
  if(face==3); La = L3; Lb = L1; Lc = L4; Ld = L2; end; 
  if(face==4); La = L4; Lb = L1; Lc = L3; Ld = L2; end;
 
  % compute warp tangential to face
  [warp1 warp2] = WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld); 
  
  blend = Lb.*Lc.*Ld;   % compute volume blending

  denom = (Lb+.5*La).*(Lc+.5*La).*(Ld+.5*La);   % modify linear blend
  ids = find(denom>tol);
  blend(ids) = (1+(alpha.*La(ids)).^2).*blend(ids)./denom(ids);

  % compute warp & blend
  shift = shift+(blend.*warp1)*t1(face,:) + (blend.*warp2)*t2(face,:);

  % fix face warp
  ids = find(La<tol & ( (Lb>tol) + (Lc>tol) + (Ld>tol) < 3));

  shift(ids,:) = warp1(ids)*t1(face,:) + warp2(ids)*t2(face,:);
end
XYZ = XYZ + shift;
X = XYZ(:,1); Y = XYZ(:,2); Z = XYZ(:,3);

[X, Y, Z] = xyztorst(X, Y, Z);
X = 0.5*(X+1);
Y = 0.5*(Y+1);
Z = 0.5*(Z+1);
plocal=[X Y Z];

function [r, s, t] = xyztorst(X, Y, Z)
    
v1 = [-1,-1/sqrt(3), -1/sqrt(6)]; v2 = [ 1,-1/sqrt(3), -1/sqrt(6)];
v3 = [ 0, 2/sqrt(3), -1/sqrt(6)]; v4 = [ 0, 0/sqrt(3),  3/sqrt(6)];

% back out right tet nodes
rhs = [X';Y';Z'] - 0.5*(v2'+v3'+v4'-v1')*ones(1,length(X));
A = [0.5*(v2-v1)',0.5*(v3-v1)',0.5*(v4-v1)'];
RST = A\rhs;
r = RST(1,:)'; s = RST(2,:)'; t = RST(3,:)';


function [X,Y,Z] = EquiNodes3D(N)

% total number of nodes
Np = (N+1)*(N+2)*(N+3)/6;

% 2) create equidistributed nodes on equilateral triangle
X = zeros(Np,1); Y = zeros(Np,1); Z = zeros(Np,1); 

sk = 1;
for n=1:N+1
  for m=1:N+2-n
    for q=1:N+3-n-m
      X(sk) = -1 + (q-1)*2/N; Y(sk) = -1 + (m-1)*2/N; Z(sk) = -1 + (n-1)*2/N;
      sk = sk+1;
    end
  end
end


function [warpx, warpy] = WarpShiftFace3D(p,pval, ~, ~,L2,L3,L4)
[dtan1,dtan2] = evalshift(p, pval, L2, L3, L4);
warpx = dtan1; warpy = dtan2;


 function [dx, dy] = evalshift(p, pval, L1, L2, L3)  
% 1) compute Gauss-Lobatto-Legendre node distribution
gaussX = -jacobigl(0,0,p);

% 3) compute blending function at each node for each edge
blend1 = L2.*L3; blend2 = L1.*L3; blend3 = L1.*L2;

% 4) amount of warp for each node, for each edge
warpfactor1 = 4*evalwarp(p, gaussX, L3-L2); 
warpfactor2 = 4*evalwarp(p, gaussX, L1-L3); 
warpfactor3 = 4*evalwarp(p, gaussX, L2-L1); 

% 5) combine blend & warp
warp1 = blend1.*warpfactor1.*(1 + (pval*L1).^2);
warp2 = blend2.*warpfactor2.*(1 + (pval*L2).^2);
warp3 = blend3.*warpfactor3.*(1 + (pval*L3).^2);

% 6) evaluate shift in equilateral triangle
dx = 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
dy = 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;


function warp = evalwarp(p, xnodes, xout)

warp = zeros(size(xout));

for i=1:p+1
  xeq(i) = -1 + 2*(p+1-i)/p;
end


for i=1:p+1
  d = (xnodes(i)-xeq(i));
  for j=2:p
    if(i~=j)
    d = d.*(xout-xeq(j))/(xeq(i)-xeq(j));
    end
  end

  if(i~=1)
    d = -d/(xeq(i)-xeq(1));
  end

  if(i~=(p+1))
    d = d/(xeq(i)-xeq(p+1));
  end

  warp = warp+d;
end

function [x] = jacobigl(alpha,beta,porder)

x = zeros(porder+1,1);
if (porder==1), x(1)=-1.0; x(2)=1.0; return; end;

xint = jacobigq(alpha+1,beta+1,porder-2);
x = [-1, xint', 1]';


function [x,w] = jacobigq(alpha,beta,N)

if (N==0), x(1)=(alpha-beta)/(alpha+beta+2); w(1) = 2; return; end;

% Form symmetric matrix from recurrence.
%J = zeros(N+1);
h1 = 2*(0:N)+alpha+beta;
J = diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1) + ...
    diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta).*...
    ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
if (alpha+beta<10*eps), J(1,1)=0.0;end;
J = J + J';

%Compute quadrature by eigenvalue solve
[V,D] = eig(J); x = diag(D);
w = (V(1,:)').^2*2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);



