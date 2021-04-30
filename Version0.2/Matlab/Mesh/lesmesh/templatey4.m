function [p,t] = templatey4(x1, y1, z1, x2, y2, z2)

xm = 0.5*(x1+x2);
ym = 0.5*(y1+y2);
zm = 0.5*(z1+z2);

xp = [x1 x2 x2 x1 x1 x2 x2 x1 xm x2 xm x1 xm x1 xm x2 x2 xm x1 xm xm];    
yp = [y1 y1 y1 y1 y2 y2 y2 y2 y1 y1 y1 y1 y1 ym ym ym ym ym ym y2 y2];        
zp = [z1 z1 z2 z2 z1 z1 z2 z2 z1 zm z2 zm zm zm zm zm z2 z2 z2 z1 z2];    

p = [xp' yp' zp'];
t = [1 9 13 12 5 20 15 14;...
     12 13 11 4 14 15 18 19;...
     14 15 18 19 5 20 21 8;...
     9 2 10 13 20 6 16 15;...
     13 10 3 11 15 16 17 18;...
     15 16 17 18 20 6 7 21];

% ptplot(p,t);