function mesh = mkmesh_sd7003(porder,elemtype,gridNum)

if nargin<1, porder=1;   end
if nargin<2, elemtype=1; end
if nargin<3, gridNum=1;  end

if gridNum==1
   n1=28*porder+1; n2=14*porder+1; n3=18*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*10);   
elseif gridNum==2
   n1=18*porder+1; n2=16*porder+1; n3=22*porder+1; 
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*10);      
elseif gridNum==3
   n1=14*porder+1; n2=6*porder+1; n3=10*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20/4, 20/4, 20/4, 20/4, 20/4, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*5);   
elseif gridNum==4
   n1=20*porder+1; n2=8*porder+1; n3=12*porder+1;    
   [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                        [20/4, 20/2, 20/2, 20/2, 20/4, 1, 1, 1, 1, 1, 1], ...
                        [10, 10, 10, 10, 10, 10, 10]*10);                       
end

[xf,yf] = read_foil('sd7003foil.txt');
[xm, ym] = cmeshmap(xf, yf, x, y, 6, 4);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<2','true'};
mesh = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);
