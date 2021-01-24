function gencode(app)

disp("generate code...");
if ~exist(char("app"), 'dir')
    mkdir(char("app"));
else
   if exist("app/opuApp.a", 'file') == 2 
       delete('app/opuApp.a');
   end
   if exist("app/cpuApp.a", 'file') == 2 
       delete('app/cpuApp.a');       
   end   
   if exist("app/gpuApp.a", 'file') == 2 
       delete('app/gpuApp.a');
   end   
end


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
    %f = pde.flux(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.flux(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Flux" + strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    error("pde.flux is not defined");
end
if isfield(pde, 'source')
    %f = pde.source(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Source" + strn, f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem("Source" + strn);
end
if isfield(pde, 'mass')
    %f = pde.mass(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Tdfunc"  + strn, f, xdg, udg, odg, wdg, uinf, param, time);
else    
    if app.model=="ModelW" || app.model == "modelW" || app.tdep==1
        error("pde.mass is not defined");
    else        
        nocodeelem("Tdfunc" + strn);
    end                
end
if isfield(pde, 'avfield')
    %f = pde.avfield(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Avfield" + strn, f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem2("Avfield" + strn);
end
if isfield(pde, 'output')
    %f = pde.output(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Output" + strn, f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem2("Output" + strn);
end
if isfield(pde, 'fbou')    
    %f = pde.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);        
    f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Fbou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("pde.fbou is not defined");
end
if isfield(pde, 'ubou')
    %f = pde.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Ubou" + strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("pde.ubou is not defined");
end
if isfield(pde, 'fhat')    
    %f = pde.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Fhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Fhat" + strn);
end
if isfield(pde, 'uhat')
    %f = pde.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Uhat" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Uhat" + strn);
end
if isfield(pde, 'stab')
    %f = pde.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Stab" + strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Stab" + strn);
end
if isfield(pde, 'initu')
    udg = pde.initu(xdg, param, uinf);
    gencodeelem3("Initu" + strn, udg, xdg, uinf, param);
else
    error("pde.initu is not defined");
end 
if isfield(pde, 'initw')
    wdg = pde.initw(xdg, param, uinf);
    gencodeelem3("Initwdg" + strn, wdg, xdg, uinf, param);
else
    if app.model=="ModelW" || app.model == "modelW"
        error("pde.initw is not defined");
    else            
        nocodeelem3("Initwdg" + strn);
    end
end
if isfield(pde, 'initv')
    odg = pde.initv(xdg, param, uinf);
    gencodeelem3("Initodg" + strn, odg, xdg, uinf, param);
else
    nocodeelem3("Initodg" + strn);
end
if isfield(pde, 'initq')
    u = pde.initu(xdg, param, uinf);
    q = pde.initq(xdg, param, uinf);
    udg = [u(:); q(:)];
    gencodeelem3("Initudg" + strn, udg, xdg, uinf, param);
    gencodeelem3("Initq" + strn, udg, xdg, uinf, param);
else
    if app.model=="ModelW" || app.model == "modelW"
        error("pde.initq is not defined");
    else
        nocodeelem3("Initudg" + strn);
        nocodeelem3("Initq" + strn);
    end
end 

end
