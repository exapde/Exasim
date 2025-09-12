function app = initializeapp(version)

app.codename = "Exasim";
app.version = version;
app.appname = "app";
app.platform = "cpu";
app.cpucompiler = "g++";
app.mpicompiler = "mpicxx";
app.gpucompiler = "nvcc";
app.mpirun = "mpirun";
%app.cpuflags = "-O2 -Wall -ldl -lm -lblas -llapack";
app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
app.gpuflags = "-lcudart -lcublas";
app.Flux = [];
app.Source = [];
app.Mass = [];
app.Avfield = [];
app.Output = [];
app.Ubou = [];
app.Fbou = [];
app.Uhat = [];
app.Fhat = [];
app.Stab = [];
app.preprocessmode = 1;
app.mpiprocs = 1;
app.nd = 1;
app.nc = 1;
app.ncu = 1;
app.ncq = 0;
app.ncp = 0;
app.nco = 0;
app.nch = 1;
app.ncx = 1;
app.ncw = 0;
app.nce = 0;
app.neb = 512*8;
app.nfb = 512*32;
app.elemtype = 1;
app.nodetype = 1;
app.pdemodel="ModelD";
app.hybrid = 0;
app.tdep = 0;
app.wave = 0;
app.linearproblem = 0;
app.debugmode = 0;
app.stgNmode = 0;
app.porder = 1;
app.pgauss = 2;
app.temporalscheme = 0;
app.torder = 1;
app.nstage = 1;
app.convStabMethod = 0;
app.diffStabMethod = 0;
app.rotatingFrame = 0;
app.viscosityModel = 0;
app.SGSmodel = 0;
app.ALE = 0;
app.AV = 0;
app.nonlinearsolver = 0;
app.linearsolver = 0;
app.NLiter = 20;
app.linearsolveriter = 200;
app.GMRESrestart = 25;
app.GMRESortho = 0;
app.preconditioner = 0;
app.precMatrixType = 0;
app.NLMatrixType = 0;
app.runmode = 0;
app.tdfunc = 1;
app.source = 1;
app.matvecorder = 1;
app.RBdim = 5;
app.saveSolFreq = 1;
app.saveSolOpt = 1;
app.timestepOffset = 0;
app.saveSolBouFreq = 0;
app.ibs = 0;
app.compudgavg = 0;
app.extFhat = 0;
app.extUhat = 0;

app.time = 0.0;
app.NLparam = 0.0;
app.NLtol = 1e-6;
app.linearsolvertol = 1e-3;
app.matvectol = 1e-6;

app.flag = [0 0];
app.problem = [0 0];
app.boundaryconditions = [0 0];
app.stgib = [0 0];

app.tau = 1.0; % stabilization parameters
app.externalparam = [0.0 0.0]; % external parameters
app.dt = [0.0 0.0];  % time steps
app.factor = [0.0 0.0];  % factors
app.physicsparam = [0.0 0.0]; % physical parameters
app.solversparam = [0.0 0.0]; % solvers parameters
app.stgdata = [0.0 0.0]; % synthetic turbulence
app.stgparam = [0.0 0.0]; % synthetic turbulence

app.soltime = 1;
app.visfilename = "dataout/output";   % filename for paraview output files
app.metis = "mpmetis";
app.gmsh = "gmsh";
app.paraview = "paraview";






