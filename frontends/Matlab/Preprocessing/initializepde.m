function pde = initializepde(version)

pde.discretization = "ldg";
pde.usecmake = 1;
pde.gencode = 1;
pde.buildexec = 0;
pde.cpucompiler = "g++";
pde.mpicompiler = "mpicxx";
pde.gpucompiler = "nvcc";
pde.usecmake = 1;
pde.mpirun = "mpirun";
pde.metis = "mpmetis";
pde.gmsh = "gmsh";
pde.paraview = "paraview";
pde.enzyme = [];
pde.codegenerator = "";

pde.codename = "Exasim";
cdir = pwd(); ii = strfind(cdir, pde.codename);
pde.exasimpath = cdir(1:(ii+5));
pde.buildpath = pde.exasimpath + "/build";
pde.backendpath = pde.exasimpath + "/backend";
pde.version = version;
pde.appname = "app";
pde.platform = "cpu";
pde.cpuflags = "-O2 -ldl -lm -lblas -llapack";
pde.gpuflags = "-lcudart -lcublas";
pde.cpuappflags = "";
pde.gpuappflags = "";
pde.cpulibflags = "";
pde.gpulibflags = "";
pde.pdemodel="ModelD";
pde.modelnumber = 0;

pde.Cxxpreprocessing = 1;
pde.preprocessmode = 1;
pde.mpiprocs = 1;
pde.nd = 1;
pde.nc = 1;
pde.ncu = 1;
pde.ncq = 0;
pde.ncp = 0;
pde.nco = 0;
pde.nch = 1;
pde.ncx = 1;
pde.ncw = 0;
pde.nce = 0;
pde.nsca = 0;
pde.nvec = 0;
pde.nten = 0;
pde.nbqoi = 0;
pde.nvqoi = 0;
pde.neb = 512*8;
pde.nfb = 512*16;
pde.elemtype = 1;
pde.nodetype = 1;
pde.hybrid = 0;
pde.tdep = 0;
pde.wave = 0;
pde.linearproblem = 0;
pde.subproblem = 0;
pde.debugmode = 0;
pde.stgNmode = 0;
pde.porder = 1;
pde.pgauss = 2;
pde.temporalscheme = 0;
pde.torder = 1;
pde.nstage = 1;
pde.convStabMethod = 0;
pde.diffStabMethod = 0;
pde.rotatingFrame = 0;
pde.viscosityModel = 0;
pde.SGSmodel = 0;
pde.ALE = 0;
pde.AV = 0;
pde.AVdistfunction = 0;
pde.AVsmoothingIter = 2;
pde.frozenAVflag = 1;
pde.nonlinearsolver = 0;
pde.linearsolver = 0;
pde.NLiter = 20;
pde.linearsolveriter = 200;
pde.GMRESrestart = 25;
pde.GMRESortho = 0;
pde.preconditioner = 1;
pde.precMatrixType = 0;
pde.ppdegree = 0;
pde.NLMatrixType = 0;
pde.runmode = 0;
pde.tdfunc = 1;
pde.sourcefunc = 1;
pde.matvecorder = 1;
pde.RBdim = 5;
pde.saveSolFreq = 1;
pde.saveSolOpt = 1;
pde.timestepOffset = 0;
pde.saveSolBouFreq = 0;
pde.ibs = 0;
pde.compudgavg = 0;
pde.extFhat = 0;
pde.extUhat = 0;
pde.extStab = 0;
pde.saveResNorm = 0;

pde.time = 0.0;
pde.NLparam = 0.0;
pde.NLtol = 1e-6;
pde.linearsolvertol = 1e-3;
pde.matvectol = 1e-3;

pde.flag = [];
pde.problem = [];
pde.boundaryconditions = [0 0];
pde.stgib = [];
pde.vindx = [];
pde.interfacefluxmap = [];
pde.avparam1 = [];
pde.avparam2 = [];

pde.tau = 1.0; 
pde.externalparam = [0.0 0.0]; 
pde.uinf = [0.0 0.0];
pde.dt = 0.0;  
pde.factor = [];  
pde.physicsparam = [0.0 0.0]; 
pde.solversparam = []; 
pde.stgdata = []; 
pde.stgparam = []; 

pde.soltime = 1;
pde.vistime = 1;
pde.visfilename = pde.buildpath + "/dataout/output";  
pde.viselem = [];

pde.dae_alpha = 1.0;
pde.dae_beta = 0.0;
pde.dae_gamma = 0.0;
pde.dae_epsilon = 0.0;
pde.dae_steps = 0;
pde.dae_dt = [];  % dual time steps

pde.fc_u = 1;
pde.fc_q = 1;
pde.denseblock=0;
pde.flux = "flux";
pde.fbou = "fbou";
pde.fhat = "fhat";
pde.source = "source";
pde.arg = {};






