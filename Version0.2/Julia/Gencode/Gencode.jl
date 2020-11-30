__precompile__()

module Gencode

using SymPy

#export syminit, gencode, compilecode
export syminit, gencode, compilecode, runcode, checkcompilers, setcompilers, tring2cmd

include("syminit.jl");
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

if isdefined(Main, Symbol("flux"))
    #f = Main.flux(xdg, udg, odg, wdg, uinf, param, time);
    f = Main.flux(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Flux", f, xdg, udg, odg, wdg, uinf, param, time);
else
    error("app.Flux is empty");
end
if isdefined(Main, Symbol("source"))
    #f = Main.source(xdg, udg, odg, wdg, uinf, param, time);
    f = Main.source(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Source", f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem("Source");
end
if isdefined(Main, Symbol("mass"))
    #f = Main.mass(xdg, udg, odg, wdg, uinf, param, time);
    f = Main.mass(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem("Tdfunc", f, xdg, udg, odg, wdg, uinf, param, time);
else
    if app.model=="ModelW" || app.model == "modelW" || app.tdep==1
        error("pde.mass is not defined");
    else
        nocodeelem("Tdfunc");
    end
end
if isdefined(Main, Symbol("avifield"))
    #f = Main.avfield(xdg, udg, odg, wdg, uinf, param, time);
    f = Main.avfield(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Avfield", f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem2("Avfield");
end
if isdefined(Main, Symbol("output"))
    #f = Main.output(xdg, udg, odg, wdg, uinf, param, time);
    f = Main.output(u, q, wdg, odg, xdg, time, param, uinf);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeelem2("Output", f, xdg, udg, odg, wdg, uinf, param, time);
else
    nocodeelem2("Output");
end
if isdefined(Main, Symbol("fbou"))
    #f = Main.fbou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = Main.fbou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    f = reshape(f,app.ncu,Int(length(f)/app.ncu));
    gencodeface("Fbou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    error("app.Fbou is empty");
end
if isdefined(Main, Symbol("ubou"))
    #f = Main.ubou(xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
    f = Main.ubou(u, q, wdg, odg, xdg, time, param, uinf, uhg, nlg, tau);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    f = reshape(f,app.ncu,Int(length(f)/app.ncu));
    gencodeface("Ubou", f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time);
else
    nocodeface("Ubou");
end
if isdefined(Main, Symbol("Fhat"))
    #f = Main.fhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = Main.fhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Fhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Fhat");
end
if isdefined(Main, Symbol("uhat"))
    #f = Main.uhat(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = Main.uhat(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Uhat", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Uhat");
end
if isdefined(Main, Symbol("stab"))
    #f = Main.stab(xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
    f = Main.stab(u1, q1, wdg1, odg1, xdg, time, param, uinf, uhg, nlg, tau, u2, q2, wdg2, odg2);
    if length(f)==1
        f = reshape([f],1,1);
    end
    f = f[:];
    gencodeface2("Stab", f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time);
else
    nocodeface2("Stab");
end
if isdefined(Main, Symbol("initu"))
    udg = Main.initu(xdg, param, uinf);
    if length(udg)==1
        udg = reshape([udg],1,1);
    end
    udg = udg[:];
    gencodeelem3("Initu", udg, xdg, uinf, param);
else
    error("initu is not defined");
end
if isdefined(Main, Symbol("initw"))
    wdg = Main.initw(xdg, param, uinf);
    if length(wdg)==1
        wdg = reshape([wdg],1,1);
    end
    wdg = wdg[:];
    gencodeelem3("Initwdg", wdg, xdg, uinf, param);
else
    nocodeelem3("Initwdg");
end
if isdefined(Main, Symbol("initv"))
    odg = Main.initv(xdg, param, uinf);
    if length(odg)==1
        odg = reshape([odg],1,1);
    end
    odg = odg[:];
    gencodeelem3("Initodg", odg, xdg, uinf, param);
else
    nocodeelem3("Initodg");
end
if isdefined(Main, Symbol("initq"))
    qdg = Main.initq(xdg, param, uinf);
    if length(qdg)==1
        qdg = reshape([qdg],1,1);
    end
    qdg = qdg[:];
    gencodeelem3("Initq", qdg, xdg, uinf, param);

    udg = Main.initu(xdg, param, uinf);
    if length(udg)==1
        udg = reshape([udg],1,1);
    end
    udg = udg[:];

    udg = [udg; qdg];
    gencodeelem3("Initudg", udg, xdg, uinf, param);
else
    if app.model == "ModelW" || app.model == "modelW"
        error("initq is not defined");
    else
        nocodeelem3("Initq");
        nocodeelem3("Initudg");
    end
end

# if isdefined(Main, Symbol("initq"))
#     qdg = Main.initq(xdg, uinf, param);
#     if length(qdg)==1
#         qdg = reshape([qdg],1,1);
#     end
#     qdg = qdg[:];
#     gencodeelem3("Initq", qdg, xdg, uinf, param);
# else
#     nocodeelem3("Initq");
# end
# if isdefined(Main, Symbol("inituq"))
#     udg = Main.inituq(xdg, uinf, param);
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
