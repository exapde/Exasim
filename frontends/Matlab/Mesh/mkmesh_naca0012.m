function [mesh,mesh1] = mkmesh_naca0012(porder,elemtype,gridNum)

if nargin<1, porder=1;   end
if nargin<2, elemtype=1; end
if nargin<3, gridNum=1;  end

if gridNum==1
   n1=28*porder+1; n2=14*porder+1; n3=18*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*5);   
elseif gridNum==2
   n1=12*porder+1; n2=10*porder+1; n3=18*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*1e4);      
elseif gridNum==3
   n1=14*porder+1; n2=6*porder+1; n3=10*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20/4, 20/4, 20/4, 20/4, 20/4, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*5);   
elseif gridNum==4
   n1=10*porder+1; n2=5*porder+1; n3=7*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]);        
elseif gridNum==5
   n1=20*porder+1; n2=40*porder+1; n3=20*porder+1;    
   %[x,y] = cmeshparam4( nxw, nfl, nfu, nr, sps, spr)
%   sps(id) : streamwise size control
%     sps(1) - ratio between the first and last elements in the wake 
%     sps(2) - ratio between leading edge and trailing edge element (lower)
%     sps(3) - ratio between leading edge and trailing edge element (upper)
%     sps(4) - ratio between the first and last elements in the wake (lower far field)
%     sps(5) - ratio between leading edge and trailing edge element (lower far field)
%     sps(6) - ratio between leading edge and trailing edge element (upper far field)
%     sps(7) - ratio between the first and last elements in the wake (upper far field)
%   sps(id) : radial size control
%     spr(1) - ratio between far-field and wake element (lower) 
%     spr(2) - ratio between far-field and trailing edge (lower)
%     spr(3) - ratio between far-field and leading edge
%     spr(4) - ratio between far-field and trailing edge (upper)
%     spr(5) - ratio between far-field and wake element (upper)  
   [x,y] = cmeshparam4(n1, n2, n2, n3, ...
                        [8/2, 4/2, 4/2, 20/2, 20/2, 20/2, 20/2, 1, 1, 1, 1], ...
                        [2, 10, 10, 10, 2, 10, 10]*3);                            
   %[x,y] = cmeshparam6( nxw, nflr, nflf, nfuf, nfur, nr, sps, spr)
%    n1=12*porder+1; n2=10*porder+1; n3=14*porder+1;    
%    [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
%                         [20/4, 20/4, 20/4, 20/4, 20/4, 1, 1, 1, 1, 1, 1], ...
%                         [10, 10, 10, 10, 10, 10, 10]*5);                         
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

figure(1); clf; plot(x, y, 'o');

[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'};
mesh = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);

mesh1 = cart2mesh(porder,x,y,[],{'true'},elemtype);

