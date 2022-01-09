function [p,t] = templatey1(x1, y1, z1, x2, y2, z2)

xm = 0.5*(x1+x2);
ym = 0.5*(y1+y2);
zm = 0.5*(z1+z2);

xp = [x1 x2 x2 x1 x1 x2 x2 x1 xm xm xm x2 x2 xm];    
yp = [y1 y1 y1 y1 y2 y2 y2 y2 y1 y1 ym ym ym ym];        
zp = [z1 z1 z2 z2 z1 z1 z2 z2 z1 z2 z1 z1 z2 z2];    

p = [xp' yp' zp'];
t = [1 9 10 4 5 11 14 8;...
     9 2 3 10 11 12 13 14;...
     11 12 13 14 5 6 7 8];

%ptplot(p,t);
