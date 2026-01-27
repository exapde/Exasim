function kkgenmodel(app)

id = app.builtinmodelID;
app.modelnumber = app.builtinmodelID;

kkdir = app.backendpath + "/Model/BuiltIn/model" + num2str(app.modelnumber);
strn = "";

if ~isfolder(kkdir)
    mkdir(kkdir);
end

text = fileread(char(app.backendpath + "/Model/BuiltIn/model.hpp")); 
fid = fopen(app.backendpath + "/Model/BuiltIn/model" + num2str(id) + "/model.hpp", 'w');  
t = string(text);          % ensure scalar string
t = strrep(t, "exasim_model_1", "exasim_model_" + num2str(id));
t = replace(t, newline, sprintf('\n'));  % normalize newlines
fwrite(fid, char(t), 'char');
fclose(fid);

text = fileread(char(app.backendpath + "/Model/BuiltIn/model.cpp")); 
fid = fopen(app.backendpath + "/Model/BuiltIn/model" + num2str(id) + "/model.cpp", 'w');  
t = string(text);          % ensure scalar string
t = strrep(t, "exasim_model_1", "exasim_model_" + num2str(id));
t = replace(t, newline, sprintf('\n'));  % normalize newlines
fwrite(fid, char(t), 'char');
fclose(fid);

editlibbuiltinmodel(id, app.backendpath + "/Model/BuiltIn/libbuiltinmodel.cpp");

disp("generate C++ code in " + kkdir);

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
udg1 = sym('udg1',[ncu+ncq 1]); % udg = (u, q)
udg2 = sym('udg2',[ncu+ncq 1]); % udg = (u, q)
odg = sym('odg',[nco 1]);
odg1 = sym('odg1',[nco 1]);
odg2 = sym('odg2',[nco 1]);
wdg = sym('wdg',[ncw 1]);
wdg1 = sym('wdg1',[ncw 1]);
wdg2 = sym('wdg2',[ncw 1]);
uhg = sym('uhg',[ncu 1]); % uhg = (uh_1,...,uh_ncu)
nlg = sym('nlg',[nd 1]); % nlg = (nl_1,...,nl_d)
tau = sym('tau',[ntau 1]); % tau = (tau_1,...,tau_ntau)
uinf = sym('uinf',[nuinf 1]); % uinf = (uinf_1,...,uinf_nuinf)
param = sym('param',[nparam 1]); % param = (param_1,...,param_nparam)        

pdemodel = str2func(app.modelfile);
pde = pdemodel();

u = udg(1:ncu);
u1 = udg1(1:ncu);
u2 = udg2(1:ncu);
if app.nc>app.ncu
    q = udg(ncu+1:end);
    q1 = udg1(ncu+1:end);
    q2 = udg2(ncu+1:end);
else
    q = [];
    q1 = [];
    q2 = [];
end

if app.hybrid == 0
  hdgnocodeelem("Flux" + strn, kkdir);  
  hdgnocodeelem("Source" + strn, kkdir);  
  hdgnocodeelem("EoS" + strn, kkdir);    
  hdgnocodeelem("Sourcew" + strn, kkdir);
  hdgnocodeelem2("Sourcewonly" + strn, kkdir);
  hdgnocodeface("Fbou" + strn, kkdir);
  hdgnocodeface2("Fbouonly" + strn, kkdir);
else
  if isfield(pde, 'flux')    
      f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);          
      hdggencodeelem("Flux" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);       
  else
      error("pde.flux is not defined");
  end
  if isfield(pde, 'source')
      f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
      hdggencodeelem("Source" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
  else        
      error("pde.source is not defined");
  end
  if isfield(pde, 'eos')
      f = pde.eos(u, q, wdg, odg, xdg, time, param, uinf);
      hdggencodeelem("EoS" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
  else    
      hdgnocodeelem("EoS" + strn, kkdir);    
  end
  if isfield(pde, 'sourcew')    
      f = pde.sourcew(u, q, wdg, odg, xdg, time, param, uinf);
      hdggencodeelem("Sourcew" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
      hdggencodeelem2("Sourcewonly" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
  else    
      hdgnocodeelem("Sourcew" + strn, kkdir);
      hdgnocodeelem2("Sourcewonly" + strn, kkdir);
  end
  if isfield(pde, 'fbouhdg')    
      f = pde.fbouhdg(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
      f = reshape(f,ncu,[]);
      hdggencodeface("Fbou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
      hdggencodeface2("Fbouonly" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
  else
      error("pde.fbouhdg is not defined");
  end  
  if isfield(pde, 'fint')    
      ncu12 = pde.ncu12;
      if ncu12==0
        ncu12 = 1;
      end              
      f = pde.fint(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
      f = reshape(f,ncu12,[]);
      hdggencodeface("Fint" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
      hdggencodeface2("Fintonly" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
  else   
      hdgnocodeface("Fint" + strn, kkdir);
      hdgnocodeface2("Fintonly" + strn, kkdir);
  end
  if isfield(pde, 'fext')    
      uext = sym('uext',[pde.ncuext 1]); 
      f = pde.fext(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, uext, tau);
      f = reshape(f,ncuext,[]);
      hdggencodefext("Fext" + strn, f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, kkdir);
      hdggencodefext2("Fextonly" + strn, f, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, kkdir);
  else   
      hdgnocodefext("Fext" + strn, kkdir);
      hdgnocodefext2("Fextonly" + strn, kkdir);
  end  
end

if isfield(pde, 'flux')    
    f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);    
    kkgencodeelem("Flux" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);       
else
    error("pde.flux is not defined");
end
if isfield(pde, 'source')
    f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("Source" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    error("pde.source is not defined");
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
if isfield(pde, 'qoivolume')
    f = pde.qoivolume(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("QoIvolume" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    kknocodeelem("QoIvolume" + strn, kkdir);
end
if isfield(pde, 'eos')
    f = pde.eos(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem2("EoS" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
    
    nf = length(f);
    nu = length(u);
    nw = length(wdg);
    
    dfdu = sym(zeros(nf,nu));
    for m = 1:nf
      for n = 1:nu
        dfdu(m,n) = diff(f(m),u(n));      
      end
    end
    kkgencodeelem2("EoSdu" + strn, dfdu, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
    
    dfdw = sym(zeros(nf,nw));
    for m = 1:nf
      for n = 1:length(wdg)
        dfdw(m,n) = diff(f(m),wdg(n));      
      end
    end
    kkgencodeelem2("EoSdw" + strn, dfdw, xdg, udg, odg, wdg, uinf, param, time, kkdir);    
else    
    kknocodeelem2("EoS" + strn, kkdir);
    kknocodeelem2("EoSdu" + strn, kkdir);
    kknocodeelem2("EoSdw" + strn, kkdir);
end
if isfield(pde, 'sourcew')    
    f = pde.sourcew(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem2("Sourcew" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
else    
    kknocodeelem2("Sourcew" + strn, kkdir);
end
if isfield(pde, 'mass')
    f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem("Tdfunc"  + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
else    
    kknocodeelem("Tdfunc" + strn, kkdir);
end
if isfield(pde, 'avfield')
    f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem2("Avfield" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
else    
    kknocodeelem2("Avfield" + strn, kkdir);
end
if isfield(pde, 'output')
    f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem2("Output" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
else    
    kknocodeelem2("Output" + strn, kkdir);
end
if isfield(pde, 'monitor')
    f = pde.monitor(u, q, wdg, odg, xdg, time, param, uinf);
    kkgencodeelem2("Monitor" + strn, f, xdg, udg, odg, wdg, uinf, param, time, kkdir);
else    
    kknocodeelem2("Monitor" + strn, kkdir);
end
if isfield(pde, 'fbou')    
    f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    kkgencodeface("Fbou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
else
    % disp("WARNING: fbou is not defined in the PDE model")
    error("pde.fbou is not defined");
end
if isfield(pde, 'ubou')
    f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    kkgencodeface("Ubou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
else
    % disp("WARNING: ubou is not defined in the PDE model")
    error("pde.ubou is not defined");
end
if isfield(pde, 'qoiboundary')    
    f = pde.qoiboundary(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    kkgencodeface("QoIboundary" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
else
    kknocodeface("QoIboundary" + strn, kkdir);
end
if isfield(pde, 'fhat')    
    f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    kkgencodeface2("Fhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, kkdir);
else
    kknocodeface2("Fhat" + strn, kkdir);
end
if isfield(pde, 'uhat')
    f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    kkgencodeface2("Uhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, kkdir);
else
    kknocodeface2("Uhat" + strn, kkdir);
end
if isfield(pde, 'stab')
    f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    kkgencodeface3("Stab" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, kkdir);
else
    kknocodeface2("Stab" + strn, kkdir);
end
if isfield(pde, 'initu')
    udg = pde.initu(xdg, param, uinf);
    kkgencodeelem3("Initu" + strn, udg, xdg, uinf, param, kkdir);
    kkgencodeelem4("Initu" + strn, udg, xdg, uinf, param, kkdir);
else
    error("pde.initu is not defined");
end 
if isfield(pde, 'initw')
    wdg = pde.initw(xdg, param, uinf);
    kkgencodeelem3("Initwdg" + strn, wdg, xdg, uinf, param, kkdir);
    kkgencodeelem4("Initwdg" + strn, wdg, xdg, uinf, param, kkdir);
else
    kknocodeelem3("Initwdg" + strn, kkdir);
    kknocodeelem4("Initwdg" + strn, kkdir);
end
if isfield(pde, 'initv')
    odg = pde.initv(xdg, param, uinf);
    kkgencodeelem3("Initodg" + strn, odg, xdg, uinf, param, kkdir);
    kkgencodeelem4("Initodg" + strn, odg, xdg, uinf, param, kkdir);
else
    kknocodeelem3("Initodg" + strn, kkdir);
    kknocodeelem4("Initodg" + strn, kkdir);
end
if isfield(pde, 'initq')
    u = pde.initu(xdg, param, uinf);
    q = pde.initq(xdg, param, uinf);
    udg = [u(:); q(:)];
    kkgencodeelem3("Initudg" + strn, udg, xdg, uinf, param, kkdir);
    kkgencodeelem3("Initq" + strn, udg, xdg, uinf, param, kkdir);
    kkgencodeelem4("Initudg" + strn, udg, xdg, uinf, param, kkdir);
    kkgencodeelem4("Initq" + strn, udg, xdg, uinf, param, kkdir);
else
    kknocodeelem3("Initudg" + strn, kkdir);
    kknocodeelem3("Initq" + strn, kkdir);
    kknocodeelem4("Initudg" + strn, kkdir);
    kknocodeelem4("Initq" + strn, kkdir);    
end 

end
