function writeapp(app,filename)

appname = 0;
app.flag = [app.tdep app.wave app.linearproblem app.debugmode app.matvecorder app.GMRESortho app.preconditioner app.precMatrixType app.NLMatrixType app.runmode app.tdfunc app.source app.modelnumber app.extFhat app.extUhat app.extStab app.subproblem app.flag];
app.problem = [app.hybrid appname app.temporalscheme app.torder app.nstage app.convStabMethod app.diffStabMethod app.rotatingFrame app.viscosityModel app.SGSmodel app.ALE app.AV app.linearsolver app.NLiter app.linearsolveriter app.GMRESrestart app.RBdim app.saveSolFreq app.saveSolOpt app.timestepOffset app.stgNmode app.saveSolBouFreq app.ibs app.dae_steps app.saveResNorm app.AVsmoothingInter app.frozenAVflag app.ppdegree app.problem];
app.factor = [app.time app.dae_alpha app.dae_beta app.dae_gamma app.dae_epsilon app.factor];
app.solversparam = [app.NLtol app.linearsolvertol app.matvectol app.NLparam app.solversparam];

ndims = zeros(40,1);
ndims[1] = app.mpiprocs;  # number of processors
ndims[2] = app.nd;
ndims[3] = 0;
ndims[4] = 0;
ndims[5] = 0;
ndims[6] = app.nc;
ndims[7] = app.ncu;
ndims[8] = app.ncq;
ndims[9] = app.ncp;
ndims[10] = app.nco;
ndims[11] = app.nch;
ndims[12] = app.ncx;
ndims[13] = app.nce;
ndims[14] = app.ncw;

nsize = zeros(20,1);
nsize[1] = length(ndims[:]);
nsize[2] = length(app.flag[:]);  # length of flag
nsize[3] = length(app.problem[:]); # length of physics
nsize[4] = length(app.externalparam[:]); # boundary data
nsize[5] = length(app.dt[:]); # number of time steps
nsize[6] = length(app.factor[:]); # length of factor
nsize[7] = length(app.physicsparam[:]); # number of physical parameters
nsize[8] = length(app.solversparam[:]); # number of solver parameters
nsize[9] = length(app.tau[:]); # number of solver parameters
nsize[10] = length(app.stgdata[:]);
nsize[11] = length(app.stgparam[:]);
nsize[12] = length(app.stgib[:]);
nsize[13] = length(app.vindx[:]);
nsize[14] = length(app.dae_dt[:]);

# app.nsize = nsize;
# app.ndims = ndims;

print("Writing app into file...\n");
fileID = open(filename,"w");

write(fileID,Float64(length(nsize[:])));
write(fileID,Float64.(nsize[:]));

if nsize[1]>0
    write(fileID,Float64.(ndims[:]));
end
if nsize[2]>0
    write(fileID,Float64.(app.flag[:]));
end
if nsize[3]>0
    write(fileID,Float64.(app.problem[:]));
end
if nsize[4]>0
    write(fileID,app.externalparam[:]);
end
if nsize[5]>0
    write(fileID,app.dt[:]);
end
if nsize[6]>0
    write(fileID,app.factor[:]);
end
if nsize[7]>0
    write(fileID,app.physicsparam[:]);
end
if nsize[8]>0
    write(fileID,app.solversparam[:]);
end
if nsize[9]>0
    write(fileID,app.tau[:]);
end
if nsize[10]>0
    write(fileID,app.stgdata[:]);
end
if nsize[11]>0
    write(fileID,app.stgparam[:]);
end
if nsize[12]>0
    write(fileID,app.stgib[:]);
end
if nsize[13]>0
    write(fileID,app.vindx[:].-1);
end
if nsize[14]>0
    write(fileID,app.dae_dt[:]);
end

if app.mutationflag == 1
    write(fileID, app.mutationopts["MixtureName"] * "X")
    write(fileID, app.mutationopts["StateModel"] * "X")
    write(fileID, app.mutationopts["ThermoDB"] * "X")
end
close(fileID);

return app;

end
