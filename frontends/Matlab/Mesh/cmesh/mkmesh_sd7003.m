gridNum=1;
porder=4;

% --- Mesh
if gridNum==1
   n1=32*porder+1; n2=14*porder+1; n3=14*porder+1; 
elseif gridNum==2
   n1=18*porder+1; n2=16*porder+1; n3=22*porder+1; 
end
[x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
                    [40, 20, 20, 20, 20, 1, 1, 1, 1, 1, 1], ...
                    [20, 10, 10, 10, 10, 10, 10]*10);

% if gridNum==1
%    n1=14*porder+1; n2=11*porder+1; n3=14*porder+1;
%    %n1=12*porder+1; n2=10*porder+1; n3=12*porder+1;
% end
% [x,y] = cmeshparam6(n1, n2, n2, n2, n2, n3, ...
%                     [10, 5, 5, 5, 5, 10, 1, 1, 1, 1, 10], ...
%                     [1, 1, 1, 1, 1, 1, 1]*50);

[xf,yf] = read_foil('sd7003foil');
[xm, ym] = cmeshmap(xf, yf, x, y, 6*2, 4*2);
% fix the wake gap
xm(1,1:n1) = xm(1,end:-1:end-(n1-1));
ym(1,1:n1) = ym(1,end:-1:end-(n1-1));

elemtype = 1;
bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<1','true'};
mesh = cart2mesh(porder,xm,ym,[],bndexpr,elemtype);

% % create DG mesh
% [p,t,p1]=cart2msh(porder,xm,ym);
%    
% bndexpr={'sqrt((p(:,1)-.5).^2+p(:,2).^2)<1','true'};
% 
% mesh = mkmesh(p,t,porder,bndexpr,0,0);
% mesh.dgnodes = p1;
% 
% [m,n]=size(xm);
% x=xm(1:porder:m,1:porder:n);
% y=ym(1:porder:m,1:porder:n);
% figure(1);clf;plot(x,y,'o');
% [m,n]=size(x);
% p = [x(:) y(:)];
% t = [1 2 m+2 m+1];
% t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
% t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
%  
  