function fb = fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time)

f = flux(xdg, udg, odg, wdg, uinf, param, time);
fb = f(1)*nlg(1)+f(2)*nlg(2)+tau(1)*(udg(1)-sym(0.0));


