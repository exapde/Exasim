function f = flux22(xdg, udg, odg, wdg, uinf, param, time)

q = udg(2:end);
f = param(1)*q;

