function gencode(app)

disp("generate code...");
codedir = app.backendpath + "/AppDriver/App";

[xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time] = syminit(app);
pdemodel = str2func(app.modelfile);
pde = pdemodel();

if app.modelnumber==0
    strn = "";
else
    strn = num2str(app.modelnumber);
end

ncu = app.ncu;
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

if isfield(pde, 'flux')    
    %f = pde.flux(xdg, udg, odg, wdg, uinf, param, time, codedir);
    f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Flux" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);   
else
    error("pde.flux is not defined");
end
if isfield(pde, 'source')
    %f = pde.source(xdg, udg, odg, wdg, uinf, param, time, codedir);
    f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Source" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);    
else    
    nocodeelem("Source" + strn, codedir);    
end
if isfield(pde, 'eos')
    f = pde.eos(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("EoS" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);
    
    nf = length(f);
    nu = length(u);
    nw = length(wdg);
    
    dfdu = sym(zeros(nf,nu));
    for m = 1:nf
      for n = 1:nu
        dfdu(m,n) = diff(f(m),u(n));      
      end
    end
    gencodeelem2("EoSdu" + strn, dfdu, xdg, udg, odg, wdg, uinf, param, time, codedir);    
    
    dfdw = sym(zeros(nf,nw));
    for m = 1:nf
      for n = 1:length(wdg)
        dfdw(m,n) = diff(f(m),wdg(n));      
      end
    end
    gencodeelem2("EoSdw" + strn, dfdw, xdg, udg, odg, wdg, uinf, param, time, codedir);    
else    
    nocodeelem2("EoS" + strn, codedir);
    nocodeelem2("EoSdu" + strn, codedir);
    nocodeelem2("EoSdw" + strn, codedir);
end
if isfield(pde, 'sourcew')    
    f = pde.sourcew(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Sourcew" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);
else    
    nocodeelem2("Sourcew" + strn, codedir);
end
if isfield(pde, 'mass')
    %f = pde.mass(xdg, udg, odg, wdg, uinf, param, time, codedir);
    f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Tdfunc"  + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);
else    
    if app.model=="ModelW" || app.model == "modelW" || app.tdep==1
        error("pde.mass is not defined");
    else        
        nocodeelem("Tdfunc" + strn, codedir);
    end                
end
if isfield(pde, 'avfield')
    %f = pde.avfield(xdg, udg, odg, wdg, uinf, param, time, codedir);
    f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Avfield" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);
else    
    nocodeelem2("Avfield" + strn, codedir);
end
if isfield(pde, 'output')
    %f = pde.output(xdg, udg, odg, wdg, uinf, param, time, codedir);
    f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Output" + strn, f, xdg, udg, odg, wdg, uinf, param, time, codedir);
else    
    nocodeelem2("Output" + strn, codedir);
end
if isfield(pde, 'fbou')    
    %f = pde.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, codedir);        
    f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Fbou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, codedir);
else
    error("pde.fbou is not defined");
end
if isfield(pde, 'ubou')
    %f = pde.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, codedir);
    f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Ubou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, codedir);
else
    error("pde.ubou is not defined");
end
if isfield(pde, 'fhat')    
    %f = pde.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
    f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Fhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
else
    nocodeface2("Fhat" + strn, codedir);
end
if isfield(pde, 'uhat')
    %f = pde.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
    f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Uhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
else
    nocodeface2("Uhat" + strn, codedir);
end
if isfield(pde, 'stab')
    %f = pde.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
    f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface3("Stab" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, codedir);
else
    nocodeface2("Stab" + strn, codedir);
end
if isfield(pde, 'initu')
    udg = pde.initu(xdg, param, uinf);
    gencodeelem3("Initu" + strn, udg, xdg, uinf, param, codedir);
else
    error("pde.initu is not defined");
end 
if isfield(pde, 'initw')
    wdg = pde.initw(xdg, param, uinf);
    gencodeelem3("Initwdg" + strn, wdg, xdg, uinf, param, codedir);
else
    if app.model=="ModelW" || app.model == "modelW"
        error("pde.initw is not defined");
    else            
        nocodeelem3("Initwdg" + strn, codedir);
    end
end
if isfield(pde, 'initv')
    odg = pde.initv(xdg, param, uinf);
    gencodeelem3("Initodg" + strn, odg, xdg, uinf, param, codedir);
else
    nocodeelem3("Initodg" + strn, codedir);
end
if isfield(pde, 'initq')
    u = pde.initu(xdg, param, uinf);
    q = pde.initq(xdg, param, uinf);
    udg = [u(:); q(:)];
    gencodeelem3("Initudg" + strn, udg, xdg, uinf, param, codedir);
    gencodeelem3("Initq" + strn, udg, xdg, uinf, param, codedir);
else
    if app.model=="ModelW" || app.model == "modelW"
        error("pde.initq is not defined");
    else
        nocodeelem3("Initudg" + strn, codedir);
        nocodeelem3("Initq" + strn, codedir);
    end
end 

end
