function hdggencode(app)

%kkdir = app.buildpath + "/model";
%kkdir = app.exasimpath + "/build/model";
kkdir = app.backendpath + "/Model";

[xdg, udg, ~, ~, wdg, ~, ~, odg, ~, ~, uhg, nlg, tau, uinf, param, time] = syminit(app);
pdemodel = str2func(app.modelfile);
pde = pdemodel();

if app.modelnumber==0
    strn = "";
else
    strn = num2str(app.modelnumber);
end

ncu = app.ncu;
u = udg(1:ncu);
if app.nc>app.ncu
    q = udg(ncu+1:end);
else
    q = [];
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
  if isfield(pde, 'fhathdg')    
      f = pde.fhathdg(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
      f = reshape(f,ncu,[]);
      hdggencodeface("Fhat" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
      hdggencodeface2("Fhatonly" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, kkdir);
  else   
      hdgnocodeface("Fhat" + strn, kkdir);
      hdgnocodeface2("Fhatonly" + strn, kkdir);      
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
      ncu12 = length(app.interfacefluxmap);
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
end

end
