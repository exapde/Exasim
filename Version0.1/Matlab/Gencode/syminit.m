function [xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time] = syminit(app)

nd = app.nd;
ncx = app.ncx;
ncu = app.ncu;
ncw = app.ncw;
ncq = app.ncq;
nco = app.nco;
ntau = length(app.tau);
nuinf = length(app.externalparam);
nparam = length(app.physicsparam);

time = sym('time');

xdg = sym('xdg',[ncx 1]); % xdg = (x_1,...,x_d)

udg = sym('udg',[ncu+ncq 1]); % udg = (u, q)
udg1 = sym('udg1',[ncu+ncq 1]); % udg = (u, q)
udg2 = sym('udg2',[ncu+ncq 1]); % udg = (u, q)

odg = sym('odg',[nco 1]);
odg1 = sym('odg1',[nco 1]);
odg2 = sym('odg2',[nco 1]);

wdg = sym('odg',[ncw 1]);
wdg1 = sym('odg1',[ncw 1]);
wdg2 = sym('odg2',[ncw 1]);

uhg = sym('uhg',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
nlg = sym('nlg',[nd 1]); % nlg = (nl_1,...,nl_d)

tau = sym('tau',[ntau 1]); % tau = (tau_1,...,tau_ntau)
uinf = sym('uinf',[nuinf 1]); % uinf = (uinf_1,...,uinf_nuinf)
param = sym('param',[nparam 1]); % param = (param_1,...,param_nparam)        

end

