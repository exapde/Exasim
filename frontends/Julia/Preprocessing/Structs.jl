__precompile__()

module Structs

export IntP, FloatP
export DmdStruct

const IntP = Int64
const FloatP = Float64

mutable struct DmdStruct
    facecon::Array{IntP,3}; # faces-to-DG nodes connectivities (used to obtain DG nodes on faces)
    eblks::Array{IntP,2};  # blocks of elements for parallel computation
    fblks::Array{IntP,2};  # blocks of faces for parallel computation
    nbsd::Array{IntP,2};   # neighboring subdomains (neighbors)
    elemsend::Array{IntP,2};  # elements to be sent to neighbors
    elemrecv::Array{IntP,2};  # elements to be received from neighbors
    elemsendpts::Array{IntP,2}; # markers for elements to be sent to neighbors
    elemrecvpts::Array{IntP,2}; # markers for elements to be received from neighbors
    elempart::Array{IntP,2};    # element partitions
    elempartpts::Array{IntP,2}; # markers for element partitions
    facepart::Array{IntP,2};    # element partitions
    facepartpts::Array{IntP,2}; # markers for element partitions
    facepartbnd::Array{IntP,2};
    t2f::Array{IntP,2};
    f::Array{IntP,2};
    bf::Array{IntP,2};
    f2t::Array{IntP,2};
    elemcon::Array{IntP,3};
    #elem2cpu::Array{IntP,2};    # element partitions
    DmdStruct() = new();
end

# mutable struct AppStruct
#     codename::String; # Exasim
#     version::String;  # Set version
#     appname::String;  # application name
#     platform::String; # CPU or GPU
#     cpucompiler::String; # Path to CPU compiler
#     mpicompiler::String; # Path to MPI compiler
#     gpucompiler::String; # Path to GPU compiler
#     mpirun::String;      # Path to MPI run command and MPI run options
#     cpuflags::String;    # options for CPU compiler
#     gpuflags::String;    # options for GGU compiler
#     # Flux::String;
#     # Source::String;
#     # Mass::String;
#     # Avfield::String;
#     # Output::String;
#     # Ubou::String;
#     # Fbou::String;
#     # Uhat::String;
#     # Fhat::String;
#     # Stab::String;
#     pdemodel::String;# used to indicate PDE model
#     pdemodelfile::String;# PDE model file name
#
#     preprocessmode::IntP; # preprocessing mode
#     mpiprocs::IntP; # number of MPI ranks
#     nd::IntP; # physical dimension
#     nc::IntP; # number of compoments of (u, q, p)
#     ncu::IntP;# number of compoments of (u)
#     ncq::IntP;# number of compoments of (q)
#     ncp::IntP;# number of compoments of (p)
#     nco::IntP;# number of compoments of (odg)
#     nch::IntP;# number of compoments of (uhat)
#     ncx::IntP;# number of compoments of (xdg)
#     ncw::IntP;# number of compoments of (wdg)
#     nce::IntP;# number of compoments of (output field)
#     neb::IntP;# number of element blocks for parallel computation
#     nfb::IntP;# number of face blocks for parallel computation
#     elemtype::IntP; # type of elements
#     nodetype::IntP; # type of nodes
#     hybrid::IntP; # discretization method
#     tdep::IntP; # flag for steady-state or time-dependent problem
#     wave::IntP; # flag for wave problem
#     linearproblem::IntP; # flag for linear problem
#     debugmode::IntP; # flag for debug mode
#     stgNmode::IntP; # number of synthetic turbulence generation modes
#
#     porder::IntP; # polymnomial degree
#     pgauss::IntP; # Gauss quadrature polynomial degree
#     temporalscheme::IntP; # temporal scheme
#     torder::IntP; # temporal order of accuracy
#     nstage::IntP; # number of RK stages
#     convStabMethod::IntP; # flag for convection stabilization method
#     diffStabMethod::IntP; # flag for diffusion stabilization method
#     rotatingFrame::IntP;  # flag for rotating frame
#     viscosityModel::IntP; # flag for viscosity model
#     SGSmodel::IntP; # flag for SGS model
#     ALE::IntP; # flag for ALE
#     AV::IntP; # flag for artificial viscosity
#     nonlinearsolver::IntP; # flag for nonlinear solver (Newton default)
#     linearsolver::IntP; # flag for linear solver (GMRES default)
#     NLiter::IntP; # maximum number of nonlinear iterations
#     linearsolveriter::IntP; # maximum number of linear iterations
#     GMRESrestart::IntP; # number of GMRES restarts
#     GMRESortho::IntP; # GMRES orthogonalization method
#     preconditioner::IntP; # flag for preconditioner method
#     precMatrixType::IntP; # flag for type of preconditioner matrix
#     NLMatrixType::IntP;
#     runmode::IntP; # flag for run mode
#     tdfunc::IntP; # flag for time-dependent function associated with time-derivative
#     source::IntP; # flag for source-term function
#
#     matvecorder::IntP; # flag for order of accuracy for matrix-vector multiplication
#     RBdim::IntP; # reduced-basis dimension
#     saveSolFreq::IntP; # flag for how often the solution saved in binary files
#     saveSolOpt::IntP; # option for how the solution be saved: 0 -> u only, 1 -> u and q
#     timestepOffset::IntP; # for restarting the simulation from the saved solution
#     saveSolBouFreq::IntP; # how often the solution be saved on a particular boundary
#     ibs::IntP; # the boundary on which the solution be saved
#     compudgavg::IntP; # flag if time-average solution is computed
#     extFhat::IntP;
#     extUhat::IntP;
#
#     time::FloatP; # starting time (usually 0, however >0 if restarting from the saved solution)
#     NLparam::FloatP;
#     NLtol::FloatP; # nonlinear solver tolerance
#     linearsolvertol::FloatP; # linear solver tolerance
#     matvectol::FloatP; # matrix-vector product tolerance
#
#     flag::Array{IntP,2};   # flag parameters
#     problem::Array{IntP,2};# problem parameters
#     boundaryconditions::Array{IntP,2};# a list of boundary condition numbers
#     stgib::Array{IntP,2};  # synthetic turbulence
#
#     uinf::Array{FloatP,2};    # boundary data
#     dt::Array{FloatP,2};      # time steps
#     factor::Array{FloatP,2};  # factors
#     physicsparam::Array{FloatP,2}; # physical parameters
#     solversparam::Array{FloatP,2}; # solvers parameters
#     tau::Array{FloatP,2}; # stabilization parameters
#     stgdata::Array{FloatP,2}; # synthetic turbulence
#     stgparam::Array{FloatP,2}; # synthetic turbulence
#
#     gmsh::String;
#     metis::String;
#     paraview::String;
#     visfilename::String;
#     visscalars;
#     visvectors;
#     viselem;
#     AppStruct() = new();
# end
#
# function InitializeAppStruct(app,version)
#     app.codename = "Exasim";
#     app.version = version;
#     app.appname = "app";
#     app.platform = "cpu";
#     app.cpucompiler = "g++";
#     app.mpicompiler = "mpicxx";
#     app.gpucompiler = "nvcc";
#     app.mpirun = "mpirun";
#     #app.cpuflags = "-O2 -Wall -ldl -lm -lblas -llapack";
#     app.cpuflags = "-O2 -ldl -lm -lblas -llapack";
#     app.gpuflags = "-lcudart -lcublas";
#     app.pdemodelfile = "";
#     # app.Source = "";
#     # app.Mass = "";
#     # app.Avfield = "";
#     # app.Output = "";
#     # app.Ubou = "";
#     # app.Fbou = "";
#     # app.Uhat = "";
#     # app.Fhat = "";
#     # app.Stab = "";
#     app.preprocessmode = 1;
#     app.mpiprocs = 1;
#     app.nd = 1;
#     app.nc = 1;
#     app.ncu = 1;
#     app.ncq = 0;
#     app.ncp = 0;
#     app.nco = 0;
#     app.nch = 1;
#     app.ncx = 1;
#     app.ncw = 0;
#     app.nce = 0;
#     app.neb = 512*8;
#     app.nfb = 512*32;
#     app.elemtype = 1;
#     app.nodetype = 1;
#     app.pdemodel="ModelD";
#     app.hybrid = 0;
#     app.tdep = 0;
#     app.wave = 0;
#     app.linearproblem = 0;
#     app.debugmode = 0;
#     app.stgNmode = 0;
#     app.porder = 1;
#     app.pgauss = 2;
#     app.temporalscheme = 0;
#     app.torder = 1;
#     app.nstage = 1;
#     app.convStabMethod = 0;
#     app.diffStabMethod = 0;
#     app.rotatingFrame = 0;
#     app.viscosityModel = 0;
#     app.SGSmodel = 0;
#     app.ALE = 0;
#     app.AV = 0;
#     app.nonlinearsolver = 0;
#     app.linearsolver = 0;
#     app.NLiter = 20;
#     app.linearsolveriter = 200;
#     app.GMRESrestart = 25;
#     app.GMRESortho = 0;
#     app.preconditioner = 0;
#     app.precMatrixType = 0;
#     app.NLMatrixType = 0;
#     app.runmode = 0;
#     app.tdfunc = 1;
#     app.source = 1;
#     app.matvecorder = 1;
#     app.RBdim = 5;
#     app.saveSolFreq = 1;
#     app.saveSolOpt = 1;
#     app.timestepOffset = 0;
#     app.saveSolBouFreq = 0;
#     app.ibs = 0;
#     app.compudgavg = 0;
#     app.extFhat = 0;
#     app.extUhat = 0;
#
#     app.time = 0.0;
#     app.NLparam = 0.0;
#     app.NLtol = 1e-6;
#     app.linearsolvertol = 1e-3;
#     app.matvectol = 1e-6;
#
#     app.flag = [0 0];
#     app.problem = [0 0];
#     app.boundaryconditions = [0 0];
#     app.stgib = [0 0];
#
#     app.tau = reshape([1.0],1,1); # stabilization parameters
#     app.uinf = [0.0 0.0]; # freestream values
#     app.dt = [0.0 0.0];  # time steps
#     app.factor = [0.0 0.0];  # factors
#     app.physicsparam = [0.0 0.0]; # physical parameters
#     app.solversparam = [0.0 0.0]; # solvers parameters
#     app.stgdata = [0.0 0.0]; # synthetic turbulence
#     app.stgparam = [0.0 0.0]; # synthetic turbulence
#
#     app.gmsh = "";
#     app.metis = "";
#     app.paraview = "";
#     app.visfilename = "";
#     app.visscalars = [];
#     app.visvectors = [];
#     app.viselem = [];
#     return app;
# end
#
# struct MasterStruct
#     nd::IntP; # physical dimension
#     porder::IntP; # polynomial degree
#     pgauss::IntP; # order of Gauss quadrature
#     elemtype::IntP; # type of elements
#     nodetype::IntP; # type of nodes
#     npe::IntP; # number of polynomials per element
#     npf::IntP; # number of polynomials per face
#     nge::IntP; # number of Gauss points per element
#     ngf::IntP; # number of Gauss points per face
#     np1d::IntP; # number of polynomials per 1D element
#     ng1d::IntP; # number of Gauss points per 1D element
#     perm::Array{IntP,2}; # to get face nodes from element nodes
#     shapegt::Array{FloatP,3}; # element shape functions at Gauss points (transpose)
#     shapegw::Array{FloatP,3}; # element shape functions at Gauss points multiplied by Gauss weights
#     shapfgt::Array{FloatP,3}; # face shape functions at Gauss points (transpose)
#     shapfgw::Array{FloatP,3}; # face shape functions at Gauss points multiplied by Gauss weights
#     shapent::Array{FloatP,3}; # element shape functions at nodes (transpose)
#     shapen::Array{FloatP,3};  # element shape functions at nodes
#     shapfnt::Array{FloatP,3}; # face shape functions at nodes (transpose)
#     shapfn::Array{FloatP,3};  # face shape functions at nodes
#     xpe::Array{FloatP,2}; # nodal points on master element
#     gpe::Array{FloatP,2}; # gauss points on master element
#     gwe::Array{FloatP,2}; # gauss weighs on master element
#     xpf::Array{FloatP,2}; # nodal points on master face
#     gpf::Array{FloatP,2}; # gauss points on master face
#     gwf::Array{FloatP,2}; # gauss weighs on master face
#     shap1dgt::Array{FloatP,3}; # 1D shape functions at Gauss points (transpose)
#     shap1dgw::Array{FloatP,3}; # 1D shape functions at Gauss points multiplied by Gauss weights
#     shap1dnt::Array{FloatP,3}; # 1D shape functions at nodes (transpose)
#     shap1dn::Array{FloatP,3}; # 1D shape functions at nodes
#     xp1d::Array{FloatP,2}; # node points on 1D element
#     gp1d::Array{FloatP,2}; # gauss points on 1D element
#     gw1d::Array{FloatP,2}; # gauss weights on 1D element
#     telem::Array{IntP,2}; # element
#     tface::Array{IntP,2}; # face
# end

# mutable struct MeshStruct
#     p::Array{FloatP,2};       # points of a linear mesh
#     t::Array{IntP,2};         # elements of a linear mesh
#     f::Array{IntP,2};         # faces of a linear mesh
#     dgnodes::Array{FloatP,3}; # spatial nodes of a high-order mesh
#     tprd::Array{IntP,2};      # elements for periodic conditions
#     boundaryexpr;             # expressions to determine boundaries
#     periodicexpr;             # expressions to determine periodic boundaries
#     curvedboundary::Array{IntP,2}; # flags for curved boundaries
#     curvedboundaryexpr;       # expressions to determine curved boundaries
#     MeshStruct() = new();
# end
#
# mutable struct SolStruct
#     UDG::Array{FloatP,3}; # intial solution (u, q) at DG nodes
#     ODG::Array{FloatP,3}; # auxilary field at DG nodes
#     WDG::Array{FloatP,3}; # dw/dt = u (wave problem) at DG nodes
#     Uout::Array{FloatP,3}; # numerical solution (u, q, w)
#     SolStruct() = new();
# end

end
