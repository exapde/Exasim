__precompile__()

module Preprocessing

using Types, Master, SymPy

export createhighordermesh, preprocessing, initializeexasim, findexec

include("xiny.jl");
include("createdgnodes.jl");
include("getelemface.jl");
include("localbasis.jl");
include("sortcolumns.jl");
include("sortrows.jl");
include("uniquerows.jl");
include("mkf2e.jl");
include("mkf2t.jl");
include("mkv2t.jl");
include("node2elem.jl");
include("mkent2elem.jl");
include("mkelconcg.jl");
include("facenumbering.jl");
include("createhighordermesh.jl");
include("partition.jl");
include("neighboringelements.jl");
include("permindex.jl");
include("elementpartition2.jl");
include("faceconnectivity2.jl");
include("facepartition2.jl");
include("meshpartition2.jl");
include("elementpartitionhdg.jl");
include("facepartitionhdg.jl");
include("meshpartitionhdg.jl");
# include("elementpartition.jl");
# include("faceconnectivity.jl");
# include("facepartition.jl");
# include("meshpartition.jl");
include("mkcgent2dgent.jl");
include("mkdge2dgf.jl");
include("mkelemblocks.jl");
include("mkfaceblocks.jl");
include("writesol.jl");
include("writeapp.jl");
include("writemaster.jl");
include("initializepde.jl");
include("initializemesh.jl");
include("initializeexasim.jl");
include("findexec.jl");

function preprocessing(app,mesh)

if app.modelnumber==0
    strn = "";
else
    strn = string(app.modelnumber);
end

model_path = joinpath(app.buildpath, "model")
datain_path = joinpath(app.buildpath, "datain", strn)
dataout_path = joinpath(app.buildpath, "dataout", strn)

# Check if directories exist and create them if they don't
if !isdir(model_path)
    mkpath(model_path)
end

if !isdir(datain_path)
    mkpath(datain_path)
end

if !isdir(dataout_path)
    mkpath(dataout_path)
end

filename = joinpath(app.buildpath, "datain", strn) * "/"
fileapp = filename * "app.bin"
filemaster = filename * "master.bin"    

if app.preprocessmode==0
    # update app structure
    app = writeapp(app,fileapp);
    return app;
end

app.nd  = size(mesh.p,1);
app.ncx = app.nd; #size(mesh.dgnodes,2);
nve,ne = size(mesh.t);

app.elemtype = 0;
if (app.nd==2) && (nve==4)
    app.elemtype=1;
end
if (app.nd==3) && (nve==8)
    app.elemtype=1;
end

if app.pgauss < 2*app.porder
  app.pgauss = 2*app.porder;
end

master = Master.mkmaster(app.nd,app.porder,app.pgauss,app.elemtype,app.nodetype);
writemaster(master,filemaster);

app.boundaryconditions = mesh.boundarycondition;
app.uinf = app.externalparam;
nuinf = length(app.uinf);
nparam = length(app.physicsparam);
xdgsym = [SymPy.symbols("xdg$i") for i=1:app.ncx];
uinfsym = [SymPy.symbols("uinf$i") for i=1:nuinf];
paramsym = [SymPy.symbols("param$i") for i=1:nparam];


if isdefined(Main, Symbol(app.modelfile))
    pdemodel = getfield(Main, Symbol(app.modelfile))
else
    pdemodel = getfield(Main, Symbol("Main"))
end
if isdefined(pdemodel, Symbol("initu"))
    udgsym = pdemodel.initu(xdgsym, paramsym, uinfsym);
    app.ncu = length(udgsym);
else
    error("pde.initu is not defined");
end
if isdefined(pdemodel, Symbol("initv"))
    vdgsym = pdemodel.initv(xdgsym, paramsym, uinfsym);
    app.nco = length(vdgsym);
else
    app.nco = 0;
end
if isdefined(pdemodel, Symbol("initw"))
    wdgsym = pdemodel.initw(xdgsym, paramsym, uinfsym);
    app.ncw = length(wdgsym);
else
    app.ncw = 0;
end

if app.model=="ModelC" || app.model=="modelC"
    app.wave = 0;
    app.nc = app.ncu;
elseif app.model=="ModelD" || app.model=="modelD"
    app.wave = 0;
    app.nc = round((app.ncu)*(app.nd+1));
elseif app.model=="ModelW" || app.model=="modelW"
    app.tdep = 1;
    app.wave = 1;
    app.nc = round((app.ncu)*(app.nd+1));
end
app.ncq = app.nc - app.ncu;
app.nch  = app.ncu;

if maximum(app.dt)>0.0
    app.tdep = 1;
else
    app.tdep = 0;
end


udgsym = [SymPy.symbols("udg$i") for i=1:app.ncu];
qdgsym = [SymPy.symbols("qdg$i") for i=1:app.ncq];
wdgsym = [SymPy.symbols("wdg$i") for i=1:app.ncw];
odgsym = [SymPy.symbols("odg$i") for i=1:app.nco];
uhatsym = [SymPy.symbols("uhg$i") for i=1:app.ncu];
nsym = [SymPy.symbols("nlg$i") for i=1:app.ncx];
time = SymPy.symbols("time");
tausym = [SymPy.symbols("tau$i") for i=1:length(app.tau)];
if isdefined(pdemodel, Symbol("visscalars")) 
    f = pdemodel.visscalars(udgsym, qdgsym, wdgsym, odgsym, xdgsym, time, paramsym, uinfsym);         
    if length(f)==1
        f = reshape([f],1,1);
    end
    app.nsca = length(f[:]);    
end
if isdefined(pdemodel, Symbol("visvectors")) 
    f = pdemodel.visvectors(udgsym, qdgsym, wdgsym, odgsym, xdgsym, time, paramsym, uinfsym);         
    if length(f)==1
        f = reshape([f],1,1);
    end
    app.nvec = length(f[:])/app.nd;    
end
if isdefined(pdemodel, Symbol("vistensors")) 
    f = pdemodel.vistensors(udgsym, qdgsym, wdgsym, odgsym, xdgsym, time, paramsym, uinfsym);         
    if length(f)==1
        f = reshape([f],1,1);
    end
    app.nten = length(f[:])/(app.nd*app.nd);    
end
if isdefined(pdemodel, Symbol("qoivolume")) 
    f = pdemodel.qoivolume(udgsym, qdgsym, wdgsym, odgsym, xdgsym, time, paramsym, uinfsym);         
    if length(f)==1
        f = reshape([f],1,1);
    end
    app.nvqoi = length(f[:]);    
end
if isdefined(pdemodel, Symbol("qoiboundary")) 
    f = pdemodel.qoiboundary(udgsym, qdgsym, wdgsym, odgsym, xdgsym, time, paramsym, uinfsym, uhatsym, nsym, tausym);             
    if length(f)==1
        f = reshape([f],1,1);
    end
    app.nbqoi = length(f[:]);    
end

print("run facenumbering...\n");
mesh.f, mesh.tprd, t2t = facenumbering(mesh.p,mesh.t,app.elemtype,mesh.boundaryexpr,mesh.periodicexpr);

mpiprocs = app.mpiprocs;
dmd = Array{DMDStruct, 1}(undef, mpiprocs);
for i = 1:mpiprocs
    dmd[i] = DMDStruct();
end

# app.neb = 50000;
# app.nfb = 100000;
#dmd = meshpartition(dmd,mesh.p,mesh.t,mesh.f,t2t,mesh.tprd,app.elemtype,app.boundaryconditions,mesh.boundaryexpr,mesh.periodicexpr,app.porder,mpiprocs,app.metis);
#dmd = meshpartition2(dmd,mesh.tprd,mesh.f,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,mpiprocs,app.metis);
if (app.hybrid ==1)
  dmd = meshpartitionhdg(dmd,mesh.tprd,mesh.f,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,mpiprocs,app.metis, app.Cxxpreprocessing);
else
  dmd = meshpartition2(dmd,mesh.tprd,mesh.f,t2t,app.boundaryconditions,app.nd,app.elemtype,app.porder,mpiprocs,app.metis, app.Cxxpreprocessing);
end

for i = 1:mpiprocs
    #xdg = mesh.dgnodes[:,:,dmd[i].elempart[:]];
    # create DG nodes
    if (size(mesh.dgnodes,3) > 0)
        xdg = mesh.dgnodes[:,:,dmd[i].elempart[:]];
    else
        xdg = createdgnodes(mesh.p,mesh.t[:,dmd[i].elempart[:]],mesh.f[:,dmd[i].elempart[:]],mesh.curvedboundary,mesh.curvedboundaryexpr,app.porder);
    end

    if mpiprocs==1
      mesh.dgnodes = xdg;        
    end
    
    if (size(mesh.udg,3) > 0)
      udg = mesh.udg[:,:,dmd[i].elempart[:]];
    else
      udg = zeros(0,0,0);
    end
    if (size(mesh.odg,3) > 0)
      odg = mesh.odg[:,:,dmd[i].elempart[:]];
    else
      odg = zeros(0,0,0);
    end
    if (size(mesh.wdg,3) > 0)
      wdg = mesh.wdg[:,:,dmd[i].elempart[:]];
    else
      wdg = zeros(0,0,0);
    end

    if mpiprocs==1
        writesol(filename  * "/sol",0,xdg, udg, odg, wdg);
    else
        writesol(filename  * "/sol",i,xdg, udg, odg, wdg);
    end

    if app.Cxxpreprocessing == 0
        cgelcon,rowent2elem,colent2elem,cgent2dgent,~ = mkcgent2dgent(xdg,1e-8);

        if mpiprocs==1
            ne = length(dmd[i].elempart);
            eblks,nbe = mkelemblocks(ne,app.neb);
            eblks = vcat(eblks, zeros(Int, 1, size(eblks,2)));
            mf = cumsum([0; dmd[i].facepartpts[:]]);
            fblks,nbf = mkfaceblocks(mf,dmd[i].facepartbnd[:],app.nfb);
            neb = maximum(eblks[2,:]-eblks[1,:])+1;
            nfb = maximum(fblks[2,:]-fblks[1,:])+1;
        else
            me = cumsum([0; dmd[i].elempartpts[1]; dmd[i].elempartpts[2]; dmd[i].elempartpts[3]]);
            eblks,nbe = mkfaceblocks(me,[0 1 2],app.neb);
            mf = cumsum([0; dmd[i].facepartpts[:]]);
            fblks,nbf = mkfaceblocks(mf,dmd[i].facepartbnd[:],app.nfb);
            neb = maximum(eblks[2,:]-eblks[1,:])+1;
            nfb = maximum(fblks[2,:]-fblks[1,:])+1;
        end

        npe = master.npe;
        nfe = size(master.perm,2);
        facecon1 = reshape(dmd[i].facecon[:,1,:],size(dmd[i].facecon,1), size(dmd[i].facecon,3));
        facecon2 = reshape(dmd[i].facecon[:,2,:],size(dmd[i].facecon,1), size(dmd[i].facecon,3));
        ind = [];
        for ii = 1:size(fblks,2)
            if fblks[3,ii]>0
                ind = vcat(ind, Vector(fblks[1,ii]:fblks[2,ii]));
            end
        end
        if length(ind)>0
            facecon2 = facecon2[:, setdiff(1:end, ind)];
        end
        rowe2f1,cole2f1,ent2ind1 = mkdge2dgf(facecon1,master.npe*length(dmd[i].elempart));
        rowe2f2,cole2f2,ent2ind2 = mkdge2dgf(facecon2,master.npe*length(dmd[i].elempart));
    end

    ndims = zeros(Int,20,1);
    ndims[1] = size(mesh.p,1);    
    ndims[2] = length(dmd[i].elempart);
    if app.Cxxpreprocessing == 0
        ndims[3] = sum(dmd[i].facepartpts);
        ndims[4] = maximum(mesh.t[:]);
        ndims[5] = nfe;
        ndims[6] = nbe;
        ndims[7] = neb;
        ndims[8] = nbf;
        ndims[9] = nfb;
    end

    nsize = zeros(Int,50,1);
    nsize[1] = length(ndims);
    if app.Cxxpreprocessing == 0
        nsize[2] = length(dmd[i].facecon);
        nsize[3] = length(eblks);
        nsize[4] = length(fblks);
    end
    nsize[5] = length(dmd[i].nbsd);
    nsize[6] = length(dmd[i].elemsend);
    nsize[7] = length(dmd[i].elemrecv);
    nsize[8] = length(dmd[i].elemsendpts);
    nsize[9] = length(dmd[i].elemrecvpts);
    nsize[10] = length(dmd[i].elempart);
    nsize[11] = length(dmd[i].elempartpts);

    if app.Cxxpreprocessing == 0    
        nsize[12] = length(cgelcon);
        nsize[13] = length(rowent2elem);
        nsize[14] = length(cgent2dgent);
        nsize[15] = length(colent2elem);
        nsize[16] = length(rowe2f1);
        nsize[17] = length(cole2f1);
        nsize[18] = length(ent2ind1);
        nsize[19] = length(rowe2f2);
        nsize[20] = length(cole2f2);
        nsize[21] = length(ent2ind2);
        nsize[22] = length(dmd[i].f2t[:])
        nsize[23] = length(dmd[i].elemcon[:])
    end
    nsize[24] = length(master.perm[:])
    nsize[25] = length(dmd[i].bf[:])
    ti = mesh.tprd[:,dmd[i].elempart[:]] .- 1;
    nsize[27] = length(ti[:]);         
    nsize[28] = length(app.boundaryconditions[:]);         

    if (mpiprocs>1)
        print("Writing mesh into file " * string(i) * "\n");
        fileID = open(string(filename  * "/mesh", string(string(i), ".bin")),"w");
    else
        print("Writing mesh into file \n");
        fileID = open(string(filename  * "/mesh", ".bin"),"w");
    end

    write(fileID,Float64(length(nsize[:])));
    write(fileID,Float64.(nsize[:]));
    write(fileID,Float64.(ndims[:]));

    if app.Cxxpreprocessing == 0
        write(fileID,Float64.(reshape(permutedims(dmd[i].facecon.-1,[2,1,3]),length(dmd[i].facecon),1)));
        write(fileID,Float64.(eblks[:]));
        write(fileID,Float64.(fblks[:]));
    end

    if (mpiprocs>1)
      write(fileID,Float64.(dmd[i].nbsd[:].-1));
      write(fileID,Float64.(dmd[i].elemsend[:].-1));
      write(fileID,Float64.(dmd[i].elemrecv[:].-1));
    end
    write(fileID,Float64.(dmd[i].elemsendpts[:]));
    write(fileID,Float64.(dmd[i].elemrecvpts[:]));
    write(fileID,Float64.(dmd[i].elempart[:].-1));
    write(fileID,Float64.(dmd[i].elempartpts[:]));

    if app.Cxxpreprocessing == 0
        write(fileID,Float64.(cgelcon[:].-1));
        write(fileID,Float64.(rowent2elem[:]));
        write(fileID,Float64.(cgent2dgent[:].-1));
        write(fileID,Float64.(colent2elem[:].-1));
        write(fileID,Float64.(rowe2f1[:]));
        write(fileID,Float64.(cole2f1[:].-1));
        write(fileID,Float64.(ent2ind1[:].-1));
        write(fileID,Float64.(rowe2f2[:]));
        write(fileID,Float64.(cole2f2[:].-1));
        write(fileID,Float64.(ent2ind2[:].-1));
        write(fileID,Float64.(dmd[i].f2t[:].-1));
        write(fileID,Float64.(dmd[i].elemcon[:].-1));
    end
    write(fileID,Float64.(master.perm[:]).-1)
    write(fileID,Float64.(dmd[i].bf[:]))
    write(fileID,Float64.(ti[:]));
    write(fileID,Float64.(app.boundaryconditions[:]));

    close(fileID);
end

app = writeapp(app,fileapp);

mesh.telem = master.telem;
mesh.tface = master.telem;
mesh.xpe = master.xpe;
mesh.xpf = master.xpf;
for i = 1:mpiprocs
    dmd[i].nbsd = reshape([0],1,1);
    #dmd[i].elem2cpu = reshape([0],1,1);
    dmd[i].elemrecv = reshape([0],1,1);
    dmd[i].elemsend = reshape([0],1,1);
    dmd[i].elemrecvpts = reshape([0],1,1);
    dmd[i].elemsendpts = reshape([0],1,1);
    dmd[i].facepartpts = reshape([0],1,1);
    dmd[i].facepartbnd = reshape([0],1,1);
    dmd[i].facecon = reshape([0],1,1,1);
end

return app, mesh, master, dmd

end

end
