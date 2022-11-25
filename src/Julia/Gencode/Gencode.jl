__precompile__()

module Gencode

using SymPy

#export syminit, gencode, compilecode
export syminit, gencode, gencodeall, compilecode, runcode, checkcompilers, setcompilers, tring2cmd

include("syminit.jl");
include("contains.jl");
include("varsassign.jl");
include("sympyassign.jl");
include("gencodebou.jl");
include("gencodeelem.jl");
include("gencodeelem2.jl");
include("gencodeelem3.jl");
include("gencodeface.jl");
include("gencodeface2.jl");
include("nocodeelem.jl");
include("nocodeelem2.jl");
include("nocodeelem3.jl");
include("nocodeface.jl");
include("nocodeface2.jl");
include("gencodeelemface.jl");
include("gencodeall.jl");
include("string2cmd.jl");
include("genlib.jl");
include("checkcompilers.jl");
include("setcompilers.jl");
include("compilecode.jl");
include("runcode.jl");

function gencode(app)

print("generate code...\n");
if !isdir("app")
    mkdir("app");
else
    if isfile("app/opuApp.a")
        rm("app/opuApp.a")
    end
    if isfile("app/cpuApp.a")
        rm("app/cpuApp.a")
    end
    if isfile("app/gpuApp.a")
        rm("app/gpuApp.a")
    end
end

xdg, udg, udg1, udg2, wdg, wdg1, wdg2, odg, odg1, odg2, uhg, nlg, tau, uinf, param, time = syminit(app);

if app.modelnumber==0
    strn = "";
else
    strn = string(app.modelnumber);
end

ncu = app.ncu;
u = udg[1:ncu];
u1 = udg1[1:ncu];
u2 = udg2[1:ncu];
if app.nc>app.ncu
    q = udg[ncu+1:end];
    q1 = udg1[ncu+1:end];
    q2 = udg2[ncu+1:end];
else
    q = [];
    q1 = [];
    q2 = [];
end

if isdefined(Main, Symbol(app.modelfile))
    pdemodel = getfield(Main, Symbol(app.modelfile))
else
    pdemodel = getfield(Main, Symbol("Main"))    
end

if isdefined(pdemodel, Symbol("flux"))
    #f = pdemodel.flux(xdg, udg, odg, wdg, uinf, param, time);
    f = pdemodel.flux(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Flux" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    error("app.Flux is empty");
end
if isdefined(pdemodel, Symbol("source"))
    #f = pdemodel.source(xdg, udg, odg, wdg, uinf, param, time);
    f = pdemodel.source(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Source" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem("Source" * strn);
end
if isdefined(pdemodel, Symbol("eos"))
    f = pdemodel.eos(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Eos" * strn, f, xdg, udg, odg, wdg, uinf, param, time);

    nf = length(f);
    nu = length(u);
    nw = length(wdg);
        
    dfdu = [SymPy.symbols("dfdu$i") for i=1:(nf*nu)];
    for n = 1:nu
      for m = 1:nf      
        dfdu[m + nf*(n-1)] = diff(f[m],u[n]);      
      end
    end        
    gencodeelem2("EoSdu" + strn, dfdu, xdg, udg, odg, wdg, uinf, param, time);    

    dfdw = [SymPy.symbols("dfdw$i") for i=1:(nf*nw)];
    for n = 1:nw
      for m = 1:nf      
        dfdw[m + nf*(n-1)] = diff(f[m],w[n]);      
      end
    end        
    gencodeelem2("EoSdw" + strn, dfdw, xdg, udg, odg, wdg, uinf, param, time);    
else
    nocodeelem2("EoS" * strn);
    nocodeelem2("EoSdu" * strn);
    nocodeelem2("EoSdw" * strn);
end
if isdefined(pdemodel, Symbol("sourcew"))
    f = pdemodel.sourcew(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Sourcew" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem2("Sourcew" * strn);
end
if isdefined(pdemodel, Symbol("mass"))
    #f = pdemodel.mass(xdg, udg, odg, wdg, uinf, param, time);
    f = pdemodel.mass(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Tdfunc" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    if app.model=="ModelW" || app.model == "modelW" || app.tdep==1
        error("pde.mass is not defined");
    else
        nocodeelem("Tdfunc" * strn);
    end
end
if isdefined(pdemodel, Symbol("avfield"))
    #f = pdemodel.avfield(xdg, udg, odg, wdg, uinf, param, time);
    f = pdemodel.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Avfield" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem2("Avfield" * strn);
end
if isdefined(pdemodel, Symbol("output"))
    #f = pdemodel.output(xdg, udg, odg, wdg, uinf, param, time);
    f = pdemodel.output(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Output" * strn, f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem2("Output" * strn);
end
if isdefined(pdemodel, Symbol("fbou"))
    #f = pdemodel.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = pdemodel.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    f = reshape(f,app.ncu,Int(length(f)/app.ncu));
    gencodeface("Fbou" * strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("app.Fbou is empty");
end
if isdefined(pdemodel, Symbol("ubou"))
    #f = pdemodel.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = pdemodel.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    f = reshape(f,app.ncu,Int(length(f)/app.ncu));
    gencodeface("Ubou" * strn, f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    nocodeface("Ubou" * strn);
end
if isdefined(pdemodel, Symbol("Fhat"))
    #f = pdemodel.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pdemodel.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Fhat" * strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Fhat" * strn);
end
if isdefined(pdemodel, Symbol("uhat"))
    #f = pdemodel.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pdemodel.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Uhat" * strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Uhat" * strn);
end
if isdefined(pdemodel, Symbol("stab"))
    #f = pdemodel.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = pdemodel.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Stab" * strn, f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Stab" * strn);
end
if isdefined(pdemodel, Symbol("initu"))
    udg = pdemodel.initu(xdg, param, uinf);
    if length(udg)==1
        udg = reshape([udg],1,1);
    end
    udg = udg[:];
    gencodeelem3("Initu" * strn, udg, xdg, uinf, param);
else
    error("initu is not defined");
end
if isdefined(pdemodel, Symbol("initw"))
    wdg = pdemodel.initw(xdg, param, uinf);
    if length(wdg)==1
        wdg = reshape([wdg],1,1);
    end
    wdg = wdg[:];
    gencodeelem3("Initwdg" * strn, wdg, xdg, uinf, param);
else
    nocodeelem3("Initwdg" * strn);
end
if isdefined(pdemodel, Symbol("initv"))
    odg = pdemodel.initv(xdg, param, uinf);
    if length(odg)==1
        odg = reshape([odg],1,1);
    end
    odg = odg[:];
    gencodeelem3("Initodg" * strn, odg, xdg, uinf, param);
else
    nocodeelem3("Initodg" * strn);
end
if isdefined(pdemodel, Symbol("initq"))
    qdg = pdemodel.initq(xdg, param, uinf);
    if length(qdg)==1
        qdg = reshape([qdg],1,1);
    end
    qdg = qdg[:];
    gencodeelem3("Initq" * strn, qdg, xdg, uinf, param);

    udg = pdemodel.initu(xdg, param, uinf);
    if length(udg)==1
        udg = reshape([udg],1,1);
    end
    udg = udg[:];

    udg = [udg; qdg];
    gencodeelem3("Initudg" * strn, udg, xdg, uinf, param);
else
    if app.model == "ModelW" || app.model == "modelW"
        error("initq is not defined");
    else
        nocodeelem3("Initq" * strn);
        nocodeelem3("Initudg" * strn);
    end
end

# if isdefined(pdemodel, Symbol("initq"))
#     qdg = pdemodel.initq(xdg, uinf, param);
#     if length(qdg)==1
#         qdg = reshape([qdg],1,1);
#     end
#     qdg = qdg[:];
#     gencodeelem3("Initq", qdg, xdg, uinf, param);
# else
#     nocodeelem3("Initq");
# end
# if isdefined(pdemodel, Symbol("inituq"))
#     udg = pdemodel.inituq(xdg, uinf, param);
#     if length(udg)==1
#         udg = reshape([udg],1,1);
#     end
#     udg = udg[:];
#     gencodeelem3("Initudg", udg, xdg, uinf, param);
# else
#     if app.model == "ModelW" || app.model == "modelW"
#         error("inituq is not defined");
#     else
#         nocodeelem3("Initudg");
#     end
# end

end

end
