function [x,y,dxds,dyds] = halfcircle(s,param)

x = cos(pi*s);
y = sin(pi*s);

dxds = -pi*sin(pi*s);
dyds =  pi*cos(pi*s);


