function syminit(app)

nd = app.nd;
nc = app.nc;
ncx = app.ncx;
ncu = app.ncu;
ncq = app.ncq;
ncp = app.ncp;
nce = app.nce;
nco = app.nco;
ncw = app.ncw;
ntau = length(app.tau);
nuinf = length(app.uinf);
nparam = length(app.physicsparam);

time = SymPy.symbols("time");

xdg = [SymPy.symbols("xdg$i") for i=1:ncx];

udg = [SymPy.symbols("udg$i") for i=1:nc];
udg1 = [SymPy.symbols("udg1$i") for i=1:nc];
udg2 = [SymPy.symbols("udg2$i") for i=1:nc];

wdg = [SymPy.symbols("wdg$i") for i=1:ncw];
wdg1 = [SymPy.symbols("wdg1$i") for i=1:ncw];
wdg2 = [SymPy.symbols("wdg2$i") for i=1:ncw];

odg = [SymPy.symbols("odg$i") for i=1:nco];
odg1 = [SymPy.symbols("odg1$i") for i=1:nco];
odg2 = [SymPy.symbols("odg2$i") for i=1:nco];

uhg = [SymPy.symbols("uhg$i") for i=1:ncu];
nlg = [SymPy.symbols("nlg$i") for i=1:nd];

tau = [SymPy.symbols("tau$i") for i=1:ntau];
uinf = [SymPy.symbols("uinf$i") for i=1:nuinf];
param = [SymPy.symbols("param$i") for i=1:nparam];

return xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time

end
