function [p,t] = refineaty(p,t,y)

% get x coordiinates
xv = unique(p(:,1));
xmin = min(xv);
xmax = max(xv); 

% get y coordiinates
yv = unique(p(:,2));
ymin = min(yv);
%ymax = max(yv); 

% find the interval [y1,y2] containing y
[~,imin] = min(abs(yv-y));
if yv(imin)>y
    i1 = imin-1;
else
    i1 = imin; 
end
y1 = yv(i1);
y2 = yv(i1+1);

% no refinement for y > y2 
[p1,t1] = removeelemement(p, t, ['y<' num2str(y2)]);
% xp = [xmin xmax xmax xmin];
% yp = [ymin ymin y2 y2];
% [p1,t1] = removeelemementinpolygon(p, t, xp, yp);
% figure(1);clf;simpplot(p1,t1);axis on;

% yes refinement for y < y1
xw = 0.5*(xv(1:end-1)+xv(2:end));
xvref = sort([xv; xw]);
[p3,t3] = quadgrid(xvref,yv(1:i1));

% transition layer between [y1, y2]
nx = length(xv)-1;
p2 = zeros(7,2,nx);
t2 = zeros(3,4,nx);
m = 1;
for i = 1:nx    
    x1 = xv(i);
    x2 = xv(i+1);
    if rem(i,2)==1
        xi = [x1 x2 x2 x1 0.5*(x1+x2) x2 0.5*(x1+x2)];    
        yi = [y1 y1 y2 y2 y1 0.5*(y1+y2) 0.5*(y1+y2)];        
        t2(:,:,i) = [1 5 7 4; 5 2 6 7; 7 6 3 4]+7*(i-1);
    else
        xi = [x1 x2 x2 x1 0.5*(x1+x2) x1 0.5*(x1+x2)];    
        yi = [y1 y1 y2 y2 y1 0.5*(y1+y2) 0.5*(y1+y2)];
        t2(:,:,i) = [1 5 7 6; 5 2 3 7; 6 7 3 4]+7*(i-1);
    end            
    p2(:,:,i) = [xi' yi'];
    m = m + 2;
end
p2 = reshape(permute(p2,[1 3 2]),[7*nx 2]);
t2 = reshape(permute(t2,[1 3 2]),[3*nx 4]);
[p2,t2] = fixmesh2(p2,t2);

% connect 3 grids 
[p12,t12] = connectmesh(p1,t1,p2,t2,1e-12);
[p,t] = connectmesh(p12,t12,p3,t3,1e-12);

% figure(1);clf;simpplot(p1,t1);axis on;
% figure(2);clf;simpplot(p2,t2);axis on;
% figure(3);clf;simpplot(p3,t3);axis on;
% pause

% [p,t] = fixmesh(p,t);



