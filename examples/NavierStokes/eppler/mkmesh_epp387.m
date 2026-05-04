
function [mesh,dgnodes] = mkmesh_epp387(porder,elemtype,gridNum)

if nargin<1, porder=1;   end
if nargin<2, elemtype=1; end
if nargin<3, gridNum=1;  end

if gridNum==-6  % p=4, very coarse mesh, 10c
   n1=24*porder+1; n2=16*porder+1; n3=20*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-8, n2-8, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30); 
elseif gridNum==-5  % p=4, Fine mesh, 50c
   n1=52*porder+1; n2=24*porder+1; n3=34*porder+1; 
   TEC = 14;
   [x,y] = cmeshparam6(n1, n2, n2-14, n2-14, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*100); 
elseif gridNum==-4  % p=2 mesh
   n1=2*38*porder+1; n2=2*24*porder+1; n3=2*32*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-12, n2-12, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30);  
elseif gridNum==-3      % p=4, Coarse mesh
   n1=30*porder+1; n2=20*porder+1; n3=26*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-10, n2-10, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30); 
elseif gridNum==-2  % p=4, Fine mesh
   n1=48*porder+1; n2=31*porder+1; n3=40*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-14, n2-14, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30); 
elseif gridNum==-1  % p=4, Medium mesh
   n1=38*porder+1; n2=24*porder+1; n3=32*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-12, n2-12, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30);   
elseif gridNum==0
   n1=46*porder+1; n2=28*porder+1; n3=40*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-14, n2-14, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30);   
elseif gridNum==1
   n1=54*porder+1; n2=32*porder+1; n3=48*porder+1; 
   TEC = 6;
   [x,y] = cmeshparam6(n1, n2, n2-10, n2-10, n2, n3, ...
                        [TEC, 1, 1, 1, 1, TEC, 1, 1, 1, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*30);   
elseif gridNum==2
   % nsp = 15 
   n1=70*porder+1; n2=48*porder+1; n3=64*porder+1; 
   TEC = 10;
   [x,y] = cmeshparam6(n1, n2, n2-40, n2-40, n2, n3, ...
                        [TEC, 1, 1/2, 1/2, 1, TEC, 1, 1/2, 1/2, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*40);                    
elseif gridNum==3
   % nsp = 20 
   n1=90*porder+1; n2=60*porder+1; n3=84*porder+1; 
   TEC = 10;
   [x,y] = cmeshparam6(n1, n2, n2-40, n2-40, n2, n3, ...
                        [TEC, 1, 1/2, 1/2, 1, TEC, 1, 1/2, 1/2, 1, TEC], ...
                        [10, 10, 10, 10, 10, 10, 10]*40);                                        
else
   error('mkmesh_epp387:UnsupportedGridNum', ...
         'Unsupported gridNum %d.', gridNum);
end

foilfile = fullfile(fileparts(mfilename('fullpath')), 'epp387_smoothed');
[xf,yf] = read_foil(foilfile);
if gridNum == -5
    Rx = 50; Ry = 50;
else
    Rx = 10; Ry = 10;
end
[xm, ym] = cmeshmap(xf, yf, x, y, Rx, Ry);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'};
mesh0 = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);

dgnodes = mesh0.dgnodes;
mesh = mesh0;
mesh.p = mesh0.p';
mesh.t = mesh0.t';
mesh.dgnodes = dgnodes;

mesh.boundaryexpr = {@(p) sqrt((p(1,:)-.5).^2+p(2,:).^2)<3, ...
                     @(p) abs(p(1,:))< 1e6 + 1e-6};
mesh.boundarycondition = [1;2];

% The C-grid mapping already supplies curved high-order DG nodes. Do not
% re-project the Eppler boundary with an approximate implicit distance.
mesh.curvedboundary = [0 0];
mesh.curvedboundaryexpr = {@(p) 0*p(1,:), @(p) 0*p(1,:)};
mesh.periodicexpr = {};
mesh.f = facenumbering(mesh.p,mesh.t,mesh.elemtype, mesh.boundaryexpr,mesh.periodicexpr);

% figure(1); clf;
% bf1 = boundaryplot(mesh,1);
% hold on;
% bf2 = boundaryplot(mesh,2);
% 
% figure(2); clf;
% hold on; 
% for i = 1:size(bf1,1)
%     p = mesh.p(:, bf1(i,:));
%     plot(p(1,:), p(2,:));
% end
% for i = 1:size(bf2,1)
%     p = mesh.p(:, bf2(i,:));
%     plot(p(1,:), p(2,:));
% end
