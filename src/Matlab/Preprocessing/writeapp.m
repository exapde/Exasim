function app = writeapp(app,filename,endian)

appname = 0;
%app.stgNmode = size(app.stgdata,1);
app.flag   = [app.tdep app.wave app.linearproblem app.debugmode app.matvecorder app.GMRESortho...  
              app.preconditioner app.precMatrixType app.NLMatrixType app.runmode app.tdfunc app.source ...
              app.modelnumber app.extFhat app.extUhat app.extStab app.subproblem app.flag];
app.problem  = [app.hybrid appname app.temporalscheme app.torder app.nstage app.convStabMethod...
               app.diffStabMethod app.rotatingFrame app.viscosityModel app.SGSmodel app.ALE app.AV...
               app.linearsolver app.NLiter app.linearsolveriter app.GMRESrestart app.RBdim ...
               app.saveSolFreq app.saveSolOpt app.timestepOffset app.stgNmode app.saveSolBouFreq app.ibs ...
               app.dae_steps app.saveResNorm app.AVsmoothingIter app.frozenAVflag app.ppdegree app.problem];
app.factor = [app.time app.dae_alpha app.dae_beta app.dae_gamma app.dae_epsilon app.factor];           
app.solversparam = [app.NLtol app.linearsolvertol app.matvectol app.NLparam app.solversparam];        

ndims = zeros(40,1);
ndims(1) = app.mpiprocs;  % number of processors
ndims(2) = app.nd;
ndims(3) = 0;
ndims(4) = 0;
ndims(5) = 0;
ndims(6) = app.nc;
ndims(7) = app.ncu;
ndims(8) = app.ncq;
ndims(9) = app.ncp;
ndims(10) = app.nco;
ndims(11) = app.nch;
ndims(12) = app.ncx;
ndims(13) = app.nce;
ndims(14) = app.ncw;
%ndims(15) = 40;
%ndims(16) = 4;

% if app.nco ~= size(app.vindx,1)
%     error("app.nco mus be equal to size(app.vindx,1)");
% end

nsize = zeros(20,1);
nsize(1) = length(ndims(:));
nsize(2) = length(app.flag(:));  % length of flag
nsize(3) = length(app.problem(:)); % length of physics
nsize(4) = length(app.externalparam(:)); 
nsize(5) = length(app.dt(:)); % number of time steps
nsize(6) = length(app.factor(:)); % length of factor
nsize(7) = length(app.physicsparam(:)); % number of physical parameters
nsize(8) = length(app.solversparam(:)); % number of solver parameters
nsize(9) = length(app.tau(:)); % number of solver parameters
nsize(10) = length(app.stgdata(:)); 
nsize(11) = length(app.stgparam(:));
nsize(12) = length(app.stgib(:));
nsize(13) = length(app.vindx(:));
nsize(14) = length(app.dae_dt(:));

app.nsize = nsize;
app.ndims = ndims;
disp('Writing app into file...'); 
fileID = fopen(filename,'w');
fwrite(fileID,length(app.nsize(:)),'double',endian);
fwrite(fileID,app.nsize(:),'double',endian);
fwrite(fileID,app.ndims(:),'double',endian);
fwrite(fileID,app.flag(:),'double',endian);
fwrite(fileID,app.problem(:),'double',endian);
fwrite(fileID,app.externalparam(:),'double',endian);
fwrite(fileID,app.dt(:),'double',endian);
fwrite(fileID,app.factor(:),'double',endian);
fwrite(fileID,app.physicsparam(:),'double',endian);
fwrite(fileID,app.solversparam(:),'double',endian);
fwrite(fileID,app.tau(:),'double',endian);
fwrite(fileID,app.stgdata(:),'double',endian);
fwrite(fileID,app.stgparam(:),'double',endian);
fwrite(fileID,app.stgib(:),'double',endian);
if (length(app.vindx(:)) > 1)    
    fwrite(fileID,app.vindx(:)-1,'double',endian);
end
fwrite(fileID,app.dae_dt(:),'double',endian);
fclose(fileID);


