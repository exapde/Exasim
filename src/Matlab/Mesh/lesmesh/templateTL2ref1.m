function [p,t] = templateTL2ref1(x0, y0, z0, x1, y1, z1)


xa = (x0 + x1)/2;
za = (z0 + z1)/2;
ya = (3*y0 + y1)/4;
yb = (y0 + y1)/2;


xp = [x0 x1 x0 x1 xa x1 x0 x0 xa x1 xa x1 xa x1 x0 xa x1 x0 xa x1 x0 xa x1];
yp = [y1 y1 y1 y1 yb yb yb yb yb yb ya ya ya ya y0 y0 y0 y0 y0 y0 y0 y0 y0];
zp = [z0 z0 z1 z1 z0 z0 za z1 z1 z1 za za z1 z1 z0 z0 z0 za za za z1 z1 z1];


p = [xp' yp' zp'];
t = [5 16 17 6 11 19 20 12;...
     5 11 12 6 9 13 14 10;...
     11 19 20 12 13 22 23 14;...
     1 15 16 5 7 18 19 11;...
     1 7 11 5 3 8 13 9;...
     7 18 19 11 8 21 22 13;...
     1 5 6 2 3 9 10 4];


%ptplot(p,t);