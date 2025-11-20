function mesh = mkmesh_naca0012(porder,elemtype,gridNum)
%[x,y] = cmeshparam6( nxw, nflr, nflf, nfuf, nfur, nr, sps, spr)
%CMESHPARAM  Creates mesh in parametric space for airfoil c-type grids
%        __________________________________ ______
%       |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%    nr |      |      |      |      |      |      |
%       |      |      |      |      |      |      |
%       |______|______|______|______|______|______|
%          nxw   nflr   nflf   nfuf   nfur
%
%   nxw  : number of subdivison in the wake
%   nflr : number of subdivision in the lower foil (rear)
%   nflf : number of subdivision in the lower foil (front)
%   nfuf : number of subdivisions in the upper foil (front)
%   nfur : number of subdivisions in the upper foil (rear)
%   nr   : number of subdivisions in the radial direction
%   sps(id) : streamwise size control
%     sps(1)  - ratio between the first and last elements in the wake 
%     sps(2)  - ratio between mid-chord and trailing edge element (lower)
%     sps(3)  - ratio between mid-chord and leading edge element (lower)
%     sps(4)  - ratio between mid-chord and leading edge element (upper)
%     sps(5)  - ratio between mid-chord and trailing edge element (upper)
%     sps(6)  - ratio between the first and last elements in the wake (lower far field)
%     sps(7)  - ratio between mid-chord and trailing edge element (lower far field)
%     sps(8)  - ratio between mid-chord and leading edge element (lower far field)
%     sps(9)  - ratio between mid-chord and leading edge element (upper far field)
%     sps(10) - ratio between mid-chord and trailing edge element (upper far field)
%     sps(11) - ratio between the first and last elements in the wake (upper far field)
%   sps(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)
%     spr(3) - ratio between far-field and mid-chord (lower)
%     spr(4) - ratio between far-field and leading edge
%     spr(5) - ratio between far-field and mid-chord (upper)
%     spr(6) - ratio between far-field and trailing edge (upper)
%     spr(7) - ratio between far-field and wake element (upper)

if nargin<1, porder=1;   end
if nargin<2, elemtype=1; end
if nargin<3, gridNum=1;  end

if gridNum==1
   n1=28*porder+1; n2=14*porder+1; n3=30*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 20, 20]*40);   
elseif gridNum==2
   %n1=14*porder+1; n2=14*porder+1; n3=20*porder+1;
   n1=22*porder+1; n2=23*porder+1; n3=30*porder+1;
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [32, 8, 2, 2, 8, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*3);      
elseif gridNum==3
   n1=2*8*porder+1; n2=2*4*porder+1; n3=2*8*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20, 10, 4, 4, 10, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*2);   
elseif gridNum==4
   n1=10*porder+1; n2=5*porder+1; n3=7*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]);                       
elseif gridNum==5
   n1=16*porder+1; n2=13*porder+1; n3=12*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [36, 10, 4, 4, 10, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*2);                       
end

thick = 12;
th = (pi:-pi/200:pi/2)';
xt = (cos(th)+1)*1.0089304129;  
xt = xt(end:-1:1);
yt=naca(xt,thick);  
xb = flipud(xt);   
yb=-naca(xb,thick);
xf =[xt; xb(2:end)];
yf =[yt; yb(2:end)];
xf(end) = xf(1);
yf(end) = yf(1);

[xm, ym] = cmeshmap(xf, yf, x, y, 8, 8);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'};
mesh = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);
mesh.p = mesh.p';
mesh.t = mesh.t';

mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, @(p) p(1,:)<7.701621+1e-6, @(p) abs(p(1,:))<20};
mesh.boundarycondition = [1;2];
mesh.curvedboundary = [1 0];
mesh.curvedboundaryexpr = {@(p) p(2,:).^2-(5*0.01*12*(0.29690*sqrt(abs(p(1,:)))-0.12600*p(1,:)-0.35160*p(1,:).^2+0.28430*p(1,:).^3-0.10150*p(1,:).^4)).^2, @(p) 0};
mesh.periodicexpr =  {};

mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype,mesh.boundaryexpr,mesh.periodicexpr);
mesh.xpe = mesh.plocal; 
mesh.telem = mesh.tlocal;

figure(1);clf;boundaryplot(mesh,1); hold on; 
boundaryplot(mesh,2); 
boundaryplot(mesh,3);

mesh.xpe = mesh.plocal;
mesh.telem = mesh.tlocal;

figure(2); clf; meshplot(mesh);
