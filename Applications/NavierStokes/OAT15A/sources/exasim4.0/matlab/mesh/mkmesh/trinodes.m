function plocal = trinodes(porder)
% TRINODES compute optimized nodal points on the master triangle
%      PLOCAL = TRINODES(PORDER)
%
%      PORDER:    Order of the complete polynomial                  
%
%      PLOCAL:    Node positions (NPL,2)
%                 NPL = (PORDER+1)*(PORDER+2)/2


alpopt = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 ...
          1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258];
          
if (porder<16)
  alpha = alpopt(porder);
else
  alpha = 5/3;
end;

Np = (porder+1)*(porder+2)/2;

% Create equidistributed nodes on equilateral triangle
L1 = zeros(Np,1); 
L3 = zeros(Np,1);
sk = 1;
for n=1:porder+1
  for m=1:porder+2-n
    L1(sk) = (n-1)/porder; L3(sk) = (m-1)/porder;
    sk = sk+1;
  end
end
L2 = 1.0-L1-L3;
x = -L2+L3; 
y = (-L2-L3+2*L1)/sqrt(3.0);

% Compute blending function at each node for each edge
blend1 = 4*L2.*L3; blend2 = 4*L1.*L3; blend3 = 4*L1.*L2;


% Amount of warp for each node, for each edge
warpf1 = warpfactor(porder,L3-L2); warpf2 = warpfactor(porder,L1-L3); warpf3 = warpfactor(porder,L2-L1);


% Combine blend & warp
warp1 = blend1.*warpf1.*(1 + (alpha*L1).^2);
warp2 = blend2.*warpf2.*(1 + (alpha*L2).^2);
warp3 = blend3.*warpf3.*(1 + (alpha*L3).^2);


% Accumulate deformations associated with each edge
x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;

L1 = (sqrt(3.0)*y+1.0)/3.0;
L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;

r = -L2 + L3 - L1; s = -L2 - L3 + L1;
x = r; 
y = s;

x = 0.5*(x+1);
y = 0.5*(y+1);
plocal = [x y];


function [x] = jacobigl(alpha,beta,porder)

x = zeros(porder+1,1);
if (porder==1), x(1)=-1.0; x(2)=1.0; return; end;

xint = jacobigq(alpha+1,beta+1,porder-2);
x = [-1, xint', 1]';

return;

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
return;

function [P] = jacobip(x,alpha,beta,porder)

% Turn points into row if needed.
xp = x;
dims = size(xp);
if (dims(2)==1), xp = xp'; end;

PL = zeros(porder+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*...
    gamma(beta+1)/gamma(alpha+beta+1);
PL(1,:) = 1.0/sqrt(gamma0);
if (porder==0), P=PL'; return; end;
gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
PL(2,:) = ((alpha+beta+2)*xp/2 + (alpha-beta)/2)/sqrt(gamma1);
if (porder==1), P=PL(porder+1,:)'; return; end;

% Repeat value in recurrence.
aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

% Forward recurrence using the symmetry of the recurrence.
for i=1:porder-1
  h1 = 2*i+alpha+beta;
  anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*...
      (i+1+beta)/(h1+1)/(h1+3));
  bnew = - (alpha^2-beta^2)/h1/(h1+2);
  PL(i+2,:) = 1/anew*( -aold*PL(i,:) + (xp-bnew).*PL(i+1,:));
  aold =anew;
end;

P = PL(porder+1,:)';
return


function [V1D] = vandermonde1d(porder,r)

% function [V1D] = Vandermonde1D(porder,r)
% Purpose : Initialize the 1D Vandermonde Matrix, V_{ij} = phi_j(r_i);

V1D = zeros(length(r),porder+1);
for j=1:porder+1
    V1D(:,j) = jacobip(r(:), 0, 0, j-1);
end;
return


function warp = warpfactor(porder, rout)

% function warp = Warpfactor(porder, rout)
% Purpose  : Compute scaled warp function at order porder based on rout interpolation nodes

% Compute LGL and equidistant node distribution
LGLr = jacobigl(0,0,porder); req  = linspace(-1,1,porder+1)';

% Compute V based on req
Veq = vandermonde1d(porder,req);

% Evaluate Lagrange polynomial at rout
Nr = length(rout); Pmat = zeros(porder+1,Nr);
for i=1:porder+1
  Pmat(i,:) = jacobip(rout, 0, 0, i-1)';
end;
Lmat = Veq'\Pmat;

% Compute warp factor
warp = Lmat'*(LGLr - req);


% Scale factor
zerof = (abs(rout)<1.0-1.0e-10); sf = 1.0 - (zerof.*rout).^2;
warp = warp./sf + warp.*(zerof-1);

return;

