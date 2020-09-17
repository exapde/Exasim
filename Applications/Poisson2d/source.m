function f = source(xdg, udg, odg, wdg, uinf, param, time)

d = 2;
x = xdg(1);
y = xdg(2);
f = (d*pi*pi)*sin(pi*x)*sin(pi*y);
