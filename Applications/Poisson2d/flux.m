function f = flux(xdg, udg, odg, wdg, uinf, param, time)

q = udg(2:end);
f = param(1)*q;

