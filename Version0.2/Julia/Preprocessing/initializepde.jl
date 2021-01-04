mutable struct PDEStruct
    codename::String; # Exasim
    version::String;  # Set version
    appname::String;  # application name
    platform::String; # CPU or GPU
    cpucompiler::String; # Path to CPU compiler
    mpicompiler::String; # Path to MPI compiler
    gpucompiler::String; # Path to GPU compiler
    mpirun::String;      # Path to MPI run command and MPI run options
    cpuflags::String;    # options for CPU compiler
    gpuflags::String;    # options for GGU compiler
    model::String;# used to indicate PDE model
    modelfile::String;# PDE model file name
    modelnumber::IntP;

    preprocessmode::IntP; # preprocessing mode
    mpiprocs::IntP; # number of MPI ranks
    nd::IntP; # physical dimension
    nc::IntP; # number of compoments of (u, q, p)
    ncu::IntP;# number of compoments of (u)
    ncq::IntP;# number of compoments of (q)
    ncp::IntP;# number of compoments of (p)
    nco::IntP;# number of compoments of (odg)
    nch::IntP;# number of compoments of (uhat)
    ncx::IntP;# number of compoments of (xdg)
    ncw::IntP;# number of compoments of (wdg)
    nce::IntP;# number of compoments of (output field)
    neb::IntP;# number of element blocks for parallel computation
    nfb::IntP;# number of face blocks for parallel computation
    elemtype::IntP; # type of elements
    nodetype::IntP; # type of nodes
    hybrid::IntP; # discretization method
    tdep::IntP; # flag for steady-state or time-dependent problem
    wave::IntP; # flag for wave problem
    linearproblem::IntP; # flag for linear problem
    debugmode::IntP; # flag for debug mode
    stgNmode::IntP; # number of synthetic turbulence generation modes

    porder::IntP; # polymnomial degree
    pgauss::IntP; # Gauss quadrature polynomial degree
    temporalscheme::IntP; # temporal scheme
    torder::IntP; # temporal order of accuracy
    nstage::IntP; # number of RK stages
    convStabMethod::IntP; # flag for convection stabilization method
    diffStabMethod::IntP; # flag for diffusion stabilization method
    rotatingFrame::IntP;  # flag for rotating frame
    viscosityModel::IntP; # flag for viscosity model
    SGSmodel::IntP; # flag for SGS model
    ALE::IntP; # flag for ALE
    AV::IntP; # flag for artificial viscosity
    nonlinearsolver::IntP; # flag for nonlinear solver (Newton default)
    linearsolver::IntP; # flag for linear solver (GMRES default)
    NLiter::IntP; # maximum number of nonlinear iterations
    linearsolveriter::IntP; # maximum number of linear iterations
    GMRESrestart::IntP; # number of GMRES restarts
    GMRESortho::IntP; # GMRES orthogonalization method
    preconditioner::IntP; # flag for preconditioner method
    precMatrixType::IntP; # flag for type of preconditioner matrix
    NLMatrixType::IntP;
    runmode::IntP; # flag for run mode
    tdfunc::IntP; # flag for time-dependent function associated with time-derivative
    source::IntP; # flag for source-term function

    matvecorder::IntP; # flag for order of accuracy for matrix-vector multiplication
    RBdim::IntP; # reduced-basis dimension
    saveSolFreq::IntP; # flag for how often the solution saved in binary files
    saveSolOpt::IntP; # option for how the solution be saved: 0 -> u only, 1 -> u and q
    timestepOffset::IntP; # for restarting the simulation from the saved solution
    saveSolBouFreq::IntP; # how often the solution be saved on a particular boundary
    ibs::IntP; # the boundary on which the solution be saved
    compudgavg::IntP; # flag if time-average solution is computed
    extFhat::IntP;
    extUhat::IntP;

    time::FloatP; # starting time (usually 0, however >0 if restarting from the saved solution)
    NLparam::FloatP;
    NLtol::FloatP; # nonlinear solver tolerance
    linearsolvertol::FloatP; # linear solver tolerance
    matvectol::FloatP; # matrix-vector product tolerance

    flag::Array{IntP,2};   # flag parameters
    problem::Array{IntP,2};# problem parameters
    boundaryconditions::Array{IntP,2};# a list of boundary condition numbers
    stgib::Array{IntP,2};  # synthetic turbulence
    vindx;  

    dt::Array{FloatP,1};      # time steps
    tau::Array{FloatP,1}; # stabilization parameters
    factor::Array{FloatP,2};  # factors
    physicsparam::Array{FloatP,2}; # physical parameters
    solversparam::Array{FloatP,2}; # solvers parameters
    stgdata::Array{FloatP,2}; # synthetic turbulence
    stgparam::Array{FloatP,2}; # synthetic turbulence
    externalparam::Array{FloatP,2};    # external parameters
    uinf::Array{FloatP,2};

    gmsh::String;
    metis::String;
    paraview::String;
    visfilename::String;
    visscalars;
    visvectors;
    viselem;
    soltime;
    visdt::FloatP;
    PDEStruct() = new();
end

function initializepde(version)
    pde = PDEStruct();

    pde.cpucompiler = "g++";
    pde.mpicompiler = "mpicxx";
    pde.gpucompiler = "nvcc";
    pde.mpirun = "mpirun";
    pde.gmsh = "gmsh";
    pde.metis = "mpmetis";
    pde.paraview = "/Applications/ParaView-5.8.1.app/Contents/MacOS/paraview";

    pde.codename = "Exasim";
    pde.version = version;
    pde.appname = "app";
    pde.platform = "cpu";
    #pde.cpuflags = "-O2 -Wall -ldl -lm -lblas -llapack";
    pde.cpuflags = "-O2 -ldl -lm -lblas -llapack";
    pde.gpuflags = "-lcudart -lcublas";
    pde.modelfile = "";
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
    pde.neb = 512*8;
    pde.nfb = 512*32;
    pde.elemtype = 1;
    pde.nodetype = 1;
    pde.model="ModelD";
    pde.modelnumber = 0;

    pde.hybrid = 0;
    pde.tdep = 0;
    pde.wave = 0;
    pde.linearproblem = 0;
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
    pde.nonlinearsolver = 0;
    pde.linearsolver = 0;
    pde.NLiter = 20;
    pde.linearsolveriter = 200;
    pde.GMRESrestart = 25;
    pde.GMRESortho = 0;
    pde.preconditioner = 0;
    pde.precMatrixType = 0;
    pde.NLMatrixType = 0;
    pde.runmode = 0;
    pde.tdfunc = 1;
    pde.source = 1;
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

    pde.time = 0.0;
    pde.NLparam = 0.0;
    pde.NLtol = 1e-6;
    pde.linearsolvertol = 1e-3;
    pde.matvectol = 1e-3;

    pde.flag = [0 0];
    pde.problem = [0 0];
    pde.boundaryconditions = [0 0];
    pde.stgib = [0 0];
    pde.vindx = [];

    pde.tau = [1.0]; # stabilization parameters
    pde.dt = [0.0];  # time steps
    pde.factor = [0.0 0.0];  # factors
    pde.physicsparam = [0.0 0.0]; # physical parameters
    pde.solversparam = [0.0 0.0]; # solvers parameters
    pde.stgdata = [0.0 0.0]; # synthetic turbulence
    pde.stgparam = [0.0 0.0]; # synthetic turbulence
    pde.externalparam = [0.0 0.0]; #
    pde.uinf = [0.0 0.0]; #

    pde.visfilename = "dataout/output";
    pde.visscalars = [];
    pde.visvectors = [];
    pde.viselem = [];
    pde.soltime = [1];
    pde.visdt = 1.0;
    return pde;
end

