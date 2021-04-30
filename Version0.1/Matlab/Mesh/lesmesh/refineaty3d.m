function [p,t] = refineaty3d(p,t,y)

% get x coordiinates
xv = unique(p(:,1));
xmin = min(xv);
xmax = max(xv); 

% get y coordiinates
yv = unique(p(:,2));
ymin = min(yv);
%ymax = max(yv); 

% get z coordiinates
zv = unique(p(:,3));
zmin = min(zv);
zmax = max(zv); 

% find the interval [y1,y2] containing y
[~,imin] = min(abs(yv-y));
if yv(imin)>y
    i1 = imin-1;
else
    i1 = imin; 
end
y1 = yv(i1-1);
y2 = yv(i1);
y3 = yv(i1+1);


% no refinement for y > y3 
[p1,t1] = removeelemement(p, t, ['y<' num2str(y3)]);

% yes refinement for y < y1
xw = 0.5*(xv(1:end-1)+xv(2:end));
xvref = sort([xv; xw]);
zw = 0.5*(zv(1:end-1)+zv(2:end));
zvref = sort([zv; zw]);
[p3,t3] = hexgrid(xvref,yv(1:i1-1),zvref);

% transition layer between [y1, y2, y3]
[p2,t2] = transitionlayer(xv, [y1, y2, y3], zv);

%p = p1; t=t1;
p = [p1; p2; p3];
t = [t1; t2+size(p1,1); t3+size(p1,1)+size(p2,1)];
[p,t]=fixmesh2(p,t);


