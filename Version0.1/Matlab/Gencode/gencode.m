function gencode(app)

disp("generate code...");
if ~exist(char("app"), 'dir')
    mkdir(char("app"));
else
   if exist("app/opuApp.a", 'file') == 2 
       delete("app/opuApp.a");
   end
   if exist("app/cpuApp.a", 'file') == 2 
       delete("app/cpuApp.a");       
   end   
   if exist("app/gpuApp.a", 'file') == 2 
       delete("app/gpuApp.a");
   end   
end

[xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time] = syminit(app);
pdemodel = str2func(app.modelfile);
pde = pdemodel();

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
    gencodeelem("Flux", f, xdg, udg, odg, wdg, uinf, param, time);
else
    error("pde.flux is not defined");
end
if isfield(pde, 'source')
    %f = pde.source(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.source(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Source", f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem("Source");
end
if isfield(pde, 'mass')
    %f = pde.mass(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.mass(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem("Tdfunc", f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem("Tdfunc");
end
if isfield(pde, 'avfield')
    %f = pde.avfield(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Avfield", f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem2("Avfield");
end
if isfield(pde, 'output')
    %f = pde.output(xdg, udg, odg, wdg, uinf, param, time);
    f = pde.output(u, q, wdg, odg, xdg, time, param, uinf);
    gencodeelem2("Output", f, xdg, udg, odg, wdg, uinf, param, time);
else    
    nocodeelem2("Output");
end
if isfield(pde, 'fbou')    
    %f = pde.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);        
    f = pde.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Fbou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("pde.fbou is not defined");
end
if isfield(pde, 'ubou')
    %f = pde.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = pde.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    f = reshape(f,ncu,[]);
    gencodeface("Ubou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("pde.ubou is not defined");
end
if isfield(pde, 'fhat')    
    %f = pde.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Fhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Fhat");
end
if isfield(pde, 'uhat')
    %f = pde.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Uhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Uhat");
end
if isfield(pde, 'stab')
    %f = pde.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pde.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    gencodeface2("Stab", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Stab");
end
if isfield(pde, 'initu')
    udg = pde.initu(xdg, param, uinf);
    gencodeelem3("Initu", udg, xdg, uinf, param);
else
    if app.model == "ModelW" || app.model == "modelW"
        nocodeelem3("Initu");        
    else        
        error("pde.initu is not defined");
    end
end 
if isfield(pde, 'initq')
    udg = pde.initq(xdg, param, uinf);
    gencodeelem3("Initq", udg, xdg, uinf, param);
else
    nocodeelem3("Initq");
end 
if isfield(pde, 'initw')
    wdg = pde.initw(xdg, param, uinf);
    gencodeelem3("Initwdg", wdg, xdg, uinf, param);
else
    nocodeelem3("Initwdg");
end
if isfield(pde, 'initv')
    odg = pde.initv(xdg, param, uinf);
    gencodeelem3("Initodg", odg, xdg, uinf, param);
else
    nocodeelem3("Initodg");
end
if isfield(pde, 'inituq')
    udg = pde.inituq(xdg, param, uinf);
    gencodeelem3("Initudg", udg, xdg, uinf, param);
else
    if app.model=="ModelW" || app.model == "modelW"
        error("pde.inituq is not defined");
    else        
        nocodeelem3("Initudg");
    end        
end 

% if ~isempty(app.Flux)
%     Flux = str2func(app.Flux);
%     f = Flux(xdg, udg, odg, wdg, uinf, param, time);    
%     gencodeelem("Flux", f, xdg, udg, odg, wdg, uinf, param, time);
% else
%     error("app.Flux is empty");
% end
% if ~isempty(app.Source)    
%     Source= str2func(app.Source);
%     f = Source(xdg, udg, odg, wdg, uinf, param, time);
%     gencodeelem("Source", f, xdg, udg, odg, wdg, uinf, param, time);
% else
%     nocodeelem("Source");
% end
% if ~isempty(app.Mass)    
%     Mass = str2func(app.Mass);
%     f = Mass(xdg, udg, odg, wdg, uinf, param, time);   
%     gencodeelem("Tdfunc", f, xdg, udg, odg, wdg, uinf, param, time);
% else
%     nocodeelem("Tdfunc");
% end
% if ~isempty(app.Avfield)
%     Avfield = str2func(app.Avfield);
%     f = Avfield(xdg, udg, odg, wdg, uinf, param, time);
%     gencodeelem2("Avfield", f, xdg, udg, odg, wdg, uinf, param, time);
% else
%     nocodeelem2("Avfield");
% end
% if ~isempty(app.Output)
%     Output = str2func(app.Output);
%     f = Output(xdg, udg, odg, wdg, uinf, param, time);
%     gencodeelem2("Output", f, xdg, udg, odg, wdg, uinf, param, time);
% else
%     nocodeelem2("Output");
% end
% if ~isempty(app.Fbou)
%     Fbou = str2func(app.Fbou);
%     f = Fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
%     gencodeface("Fbou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
% else
%     error("app.Fbou is empty");
% end
% if ~isempty(app.Ubou)
%     Ubou = str2func(app.Ubou);
%     f = Ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
%     gencodeface("Ubou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
% else
%     error("app.Ubou is empty");
% end
% if ~isempty(app.Fhat)
%     Fhat = str2func(app.Fhat);
%     f = Fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
%     gencodeface2("Fhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
% else
%     nocodeface2("Fhat");
% end
% if ~isempty(app.Uhat)
%     Uhat = str2func(app.Uhat);
%     f = Uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
%     gencodeface2("Uhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
% else
%     nocodeface2("Uhat");
% end
% if ~isempty(app.Stab)    
%     Stab = str2func(app.Stab);
%     f = Stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
%     gencodeface2("Stab", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
% else
%     nocodeface2("Stab");
% end

end
