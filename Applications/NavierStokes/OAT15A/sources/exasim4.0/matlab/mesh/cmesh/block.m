function X = block( n, sp, x)
%BLOCK  Create block coordinates.
%
%    Syntax: X = block( n, sp, x)
%
%    n(1:ndim) number of points in each direction
%      The number of components in n determines the dimension
%    sp(ned) desired sapcing distribution over each edge
%      For dim=1 -> ned=1
%          dim=2 -> ned=4
%          dim=3 -> ned=12
%      sp(i) > 0 is the ratio of the last to the first subdivision on edge i
%      sp(i) < 0 is the negative of the ratio betweend the subdivision in
%                the middle and the end subdivision (symmetric)
%    x(1:2^ndim,1:ndim) gives new coordinates of the unit cube for mapping
%
%    X(ndim,n(1),n(2),n(3)) block coordinates
%
%    EXAMPLES
%
%    Boundary layer in 2D:
%      X = block([20,20], [1,1,20,20],[0,0;3,0;0,1;3,1]);
%      surf(squeeze(X(1,:,:)),squeeze(X(2,:,:)),0*squeeze(X(1,:,:)))
%      view(2),axis equal

num = [1,4,12];
dim = prod(size(n));
if nargin < 2, sp=ones(num(dim),1); end

switch dim
case 1
    X = block1D(n,sp);
case 2
    X = block2D(n(1),n(2),sp(1:2),sp(3:4));
case 3
    X = block3D(n(1),n(2),n(3),sp(1:4),sp(5:8),sp(9:12));
end

if nargin > 2
    X = mapp(reshape(X,dim,[])',x)';
    switch dim
        case 1
            X = reshape(X,dim,n(1));
        case 2
            X = reshape(X,dim,n(1),n(2));
        case 3
            X = reshape(X,dim,n(1),n(2),n(3));
    end
end


function X = block1D( nx, spx)

if nargin<2, spx=1; end

X = zeros(1,nx);

if spx == 1 || spx == -1
   dx = 1/(nx-1);
   for j=2:nx
       X(j) = X(j-1) + dx;
   end
elseif spx > 0
   rx = spx^(1/(nx-2));
   dx = (rx-1)/(rx^(nx-1)-1);
   for j=2:nx
       X(j) = X(j-1) + dx;
       dx = dx*rx;
   end
elseif spx < 0
   nxp = 2*(nx-1)+1;
   XP = zeros(1,nxp);
   nx2 = (nxp-1)/2+1;
   rx = (-spx)^(1/(nx2-2));
   dx = 0.5*(rx-1)/(rx^(nx2-1)-1);
   for j=2:nx2
       XP(j) = XP(j-1) + dx;
       dx = dx*rx;
   end
   for j=nx2+1:nxp
       dx = dx/rx;
       XP(j) = XP(j-1) + dx;
   end
   X = XP(1:2:end);
end

function X = block2D( nx, ny, spx, spy)

X = zeros(2,nx,ny);

X(1,nx,:) = 1;
X(2,:,ny) = 1;

for i=1:2
    X(1,:,1+(i-1)*(ny-1)) = block1D( nx, spx(i));
    X(2,1+(i-1)*(nx-1),:) = block1D( ny, spy(i));
end

X = fill2D(X);


function X = block3D( nx, ny, nz, spx, spy, spz)

X = zeros(3,nx,ny,nz);

X(1,nx,:,:) = 1;
X(2,:,ny,:) = 1;
X(3,:,:,nz) = 1;

for j=1:2
  for i=1:2
      X(1,:,1+(i-1)*(ny-1),1+(j-1)*(nz-1)) = block1D( nx, spx(i+2*(j-1)));
      X(2,1+(i-1)*(nx-1),:,1+(j-1)*(nz-1)) = block1D( ny, spy(i+2*(j-1)));
      X(3,1+(i-1)*(nx-1),1+(j-1)*(ny-1),:) = block1D( nz, spz(i+2*(j-1)));
  end
end

X2 = squeeze(X([1,2],:,:,1)); X2 = fill2D(X2); X([1,2],:,:,1) = X2;
X2 = squeeze(X([1,2],:,:,nz)); X2 = fill2D(X2); X([1,2],:,:,nz) = X2;

X2 = squeeze(X([1,3],:,1,:)); X2 = fill2D(X2); X([1,3],:,1,:) = X2;
X2 = squeeze(X([1,3],:,ny,:)); X2 = fill2D(X2); X([1,3],:,ny,:) = X2;

X2 = squeeze(X([2,3],1,:,:)); X2 = fill2D(X2); X([2,3],1,:,:) = X2;
X2 = squeeze(X([2,3],nx,:,:)); X2 = fill2D(X2); X([2,3],nx,:,:) = X2;

X = fill3D(X);


function X = fill2D(X)

nx = size(X,2); ny = size(X,3);

for i=2:nx-1
    for j=2:ny-1
        X0 = X(1,i,1); Y0 = X(2,i,1);  X1 = X(1,i,ny); Y1 = X(2,i,ny);
        X2 = X(1,1,j); Y2 = X(2,1,j);  X3 = X(1,nx,j); Y3 = X(2,nx,j);
        A = [X1-X0, X2-X3; Y1-Y0, Y2-Y3]; B=[X2-X0;Y2-Y0];
        P = inv(A)*B;
        X(1,i,j) = X0 + P(1)*(X1-X0);
        X(2,i,j) = Y0 + P(1)*(Y1-Y0);
    end
end


function X = fill3D(X)

nx = size(X,2); ny = size(X,3); nz = size(X,4);

X0 = zeros(3,1); X1 = zeros(3,1); X2 = zeros(3,1); X3 = zeros(3,1);
Y0 = zeros(3,1); Y1 = zeros(3,1); Y2 = zeros(3,1); Y3 = zeros(3,1);
Z0 = zeros(3,1); Z1 = zeros(3,1); Z2 = zeros(3,1); Z3 = zeros(3,1);

for i=2:nx-1
    for j=2:ny-1
        for k=2:nz-1
            X0=squeeze(X(:,i,1,1)); X1=squeeze(X(:,i,ny,1)); X2=squeeze(X(:,i,1,nz)); X3=squeeze(X(:,i,ny,nz));      
            Y0=squeeze(X(:,1,j,1)); Y1=squeeze(X(:,nx,j,1)); Y2=squeeze(X(:,1,j,nz)); Y3=squeeze(X(:,nx,j,nz));
            Z0=squeeze(X(:,1,1,k)); Z1=squeeze(X(:,1,ny,k)); Z2=squeeze(X(:,nx,1,k)); Z3=squeeze(X(:,nx,ny,k));
            X(:,i,j,k) = fpoint3D( X0, X1, X2, X3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3);
        end
    end
end


function PX = fpoint3D( X0, X1, X2, X3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3)

P = 0.5*ones(6,1);
DP = ones(6,1);

while norm(DP) > 1.e-9
    
    PX = X0*(1-P(1))*(1-P(2)) + X1*P(1)*(1-P(2)) + X2*(1-P(1))*P(2) + X3*P(1)*P(2);
    PY = Y0*(1-P(3))*(1-P(4)) + Y1*P(3)*(1-P(4)) + Y2*(1-P(3))*P(4) + Y3*P(3)*P(4);
    PZ = Z0*(1-P(5))*(1-P(6)) + Z1*P(5)*(1-P(6)) + Z2*(1-P(5))*P(6) + Z3*P(5)*P(6);

    PXr = -X0*(1-P(2)) + X1*(1-P(2)) - X2*P(2) + X3*P(2);
    PXs = -X0*(1-P(1)) - X1*P(1) + X2*(1-P(1)) + X3*P(1);
    PYr = -Y0*(1-P(4)) + Y1*(1-P(4)) - Y2*P(4) + Y3*P(4);
    PYs = -Y0*(1-P(3)) - Y1*P(3) + Y2*(1-P(3)) + Y3*P(3);
    PZr = -Z0*(1-P(6)) + Z1*(1-P(6)) - Z2*P(6) + Z3*P(6);
    PZs = -Z0*(1-P(5)) - Z1*P(5) + Z2*(1-P(5)) + Z3*P(5);

    A = zeros(6,6);
    R = zeros(6,1);

    A(1:3,1) = PXr;
    A(4:6,1) = PXr;
    A(1:3,2) = PXs;
    A(4:6,2) = PXs;
    A(1:3,3) = -PYr;
    A(1:3,4) = -PYs;
    A(4:6,5) = -PZr;
    A(4:6,6) = -PZs;

    R(1:3) = -(PX-PY);
    R(4:6) = -(PX-PZ);

    DP = inv(A)*R;

    P = P + DP;
end
