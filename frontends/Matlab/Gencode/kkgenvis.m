function kkgenvis(app)

kkdir = app.backendpath + "/Model";
strn = "";

if ~isfolder(kkdir)
    mkdir(kkdir);
end

app.nco = app.ncv;
app.ncx = app.nd;
nd = app.nd;
ncx = app.ncx;
ncu = app.ncu;
ncw = app.ncw;
ncq = app.ncq;
nco = app.nco;
ntau = app.ntau;
nuinf = app.neta;
nparam = app.nmu;
app.nc = ncu + ncq;

time = sym('time');
xdg = sym('xdg',[ncx 1]); % xdg = (x_1,...,x_d)
udg = sym('udg',[ncu+ncq 1]); % udg = (u, q)
odg = sym('odg',[nco 1]);
wdg = sym('wdg',[ncw 1]);
uinf = sym('uinf',[nuinf 1]); % uinf = (uinf_1,...,uinf_nuinf)
param = sym('param',[nparam 1]); % param = (param_1,...,param_nparam)        

pdemodel = str2func(app.modelfile);
pde = pdemodel();

u = udg(1:ncu);
if app.nc>app.ncu
    q = udg(ncu+1:end);
else
    q = [];
end

if isfield(pde, 'visscalars')
    f = pde.visscalars(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("VisScalars" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    kknocodeelem("VisScalars" + strn, kkdir);
end
if isfield(pde, 'visvectors')
    f = pde.visvectors(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("VisVectors" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    kknocodeelem("VisVectors" + strn, kkdir);
end
if isfield(pde, 'vistensors')
    f = pde.vistensors(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("VisTensors" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    kknocodeelem("VisTensors" + strn, kkdir);
end

end
