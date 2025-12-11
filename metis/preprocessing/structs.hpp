#ifndef __STRUCTS_HPP__
#define __STRUCTS_HPP__

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct FunctionDef {
    std::string name;
    std::string output;
    int outputsize = 0;
    std::vector<std::string> args;
    std::vector<std::string> body;    
    std::unordered_map<std::string, std::pair<int, int>> matrices; // name -> (rows, cols)
};

struct ParsedSpec {
    std::vector<std::string> scalars;
    std::unordered_map<std::string, int> vectors;
    std::vector<std::string> namevectors;
    std::vector<std::string> jacobian;
    std::vector<std::string> hessian;
    std::vector<std::string> batch;
    std::vector<std::string> outputs;
    std::vector<std::string> exasimfunctions = {
        "Flux", "Source", "Tdfunc", "Ubou", "Fbou", "FbouHdg",
        "Sourcew", "Output", "Monitor", "Initu", "Initq", "Inituq",
        "Initw", "Initv", "Avfield", "Fint", "EoS", "VisScalars", 
        "VisVectors", "VisTensors", "QoIvolume", "QoIboundary"};
    std::vector<bool> isoutput;     
    std::string datatype = "dstype";
    std::string framework = "kokkos";
    std::string codeformat = "exasim";    
    std::string exasimpath = "";
    std::string modelpath = "";
    std::string symenginepath = "";
    std::string modelfile = "";
    bool exasim;
    std::vector<FunctionDef> functions;
};

// Struct to hold all parsed input parameters
struct InputParams {
    std::string pdeappfile;
    std::unordered_map<std::string, std::string> stringParams;
    std::unordered_map<std::string, double> doubleParams;
    std::unordered_map<std::string, int> intParams;

    std::vector<int> interfaceConditions;
    std::vector<int> interfaceFluxmap;
    std::vector<int> boundaryConditions;
    std::vector<int> curvedBoundaries;
    std::vector<int> periodicBoundaries1;
    std::vector<int> periodicBoundaries2;
    std::vector<int> cartGridPart;
    
    std::vector<double> dae_dt;
    std::vector<double> dt;
    std::vector<double> tau;
    std::vector<double> physicsParam;
    std::vector<double> externalParam;
    std::vector<double> vindx;
    std::vector<double> avparam1, avparam2;       
    std::vector<double> stgib;
    std::vector<double> stgdata;
    std::vector<double> stgparam;
    
    std::vector<std::string> boundaryExprs;    
    std::vector<std::string> curvedBoundaryExprs;    
    std::vector<std::string> periodicExprs1;    
    std::vector<std::string> periodicExprs2;        
    
    std::unordered_set<std::string> foundKeys;
};

struct PDE {
    std::string discretization = "ldg";
    std::string exasimpath = "";
    std::string datapath = "";
    std::string datainpath = "";
    std::string dataoutpath = "";
    std::string platform = "cpu";   
    std::string model = "ModelD";
    std::string pdeappfile = "";
    std::string modelfile = "pdemodel.txt";
    std::string meshfile = "mesh.bin";
    std::string xdgfile = "";
    std::string udgfile = "";
    std::string vdgfile = "";
    std::string wdgfile = "";
    std::string uhatfile = "";
    std::string partitionfile = "";

    int gendatain = 1;
    int gencode = 1; // 1 for code generation, 0 for no code generation
    int writemeshsol = 1; // 1 for writing mesh solution, 0 for no writing
    int modelnumber = 0;
    int mpiprocs = 1;
    int nd = 1, nc = 1, ncu = 1, ncq = 0, ncp = 0, ncv = 0;
    int nch = 1, ncx = 1, ncw = 0, nce = 0, np=0, nve=0, ne=0;
    int nsca=0, nvec=0, nten=0, nsurf=0, nvqoi=0;
    int neb = 512 * 8;
    int nfb = 512 * 16;
    int elemtype = 1;
    int nodetype = 1;
    int hybrid = 0;
    int tdep = 0;
    int wave = 0;
    int linearproblem = 0;
    int subproblem = 0;
    int debugmode = 0;
    int stgNmode = 0;
    int porder = 1;
    int pgauss = 2;
    int temporalscheme = 0;
    int torder = 1;
    int nstage = 1;
    int convStabMethod = 0;
    int diffStabMethod = 0;
    int rotatingFrame = 0;
    int viscosityModel = 0;
    int SGSmodel = 0;
    int ALE = 0;
    int AV = 0;
    int AVdistfunction = 0;
    int AVsmoothingIter = 2;
    int frozenAVflag = 1;
    int nonlinearsolver = 0;
    int linearsolver = 0;
    int NewtonIter = 20;
    int GMRESiter = 200;
    int GMRESrestart = 25;
    int GMRESortho = 0;
    int preconditioner = 0;
    int precMatrixType = 0;
    int ppdegree = 0;
    int NLMatrixType = 0;
    int runmode = 0;
    int tdfunc = 1;
    int sourcefunc = 1;
    int matvecorder = 1;
    int RBdim = 5;
    int saveSolFreq = 1;
    int saveSolOpt = 1;
    int timestepOffset = 0;
    int saveSolBouFreq = 0;
    int ibs = 0;
    int compudgavg = 0;
    int extFhat = 0;
    int extUhat = 0;
    int extStab = 0;
    int saveResNorm = 0;
    int dae_steps = 0;

    int coupledinterface = 0; 
    int coupledcondition = 0;
    int coupledboundarycondition = 0;    
    
    double time = 0.0;
    double NLparam = 0.0;
    double NewtonTol = 1e-6;
    double GMREStol = 1e-3;
    double matvectol = 1e-3;
    double dae_alpha = 1.0;
    double dae_beta = 0.0;
    double dae_gamma = 0.0;
    double dae_epsilon = 0.0;    
            
    std::vector<int> interfaceFluxmap;    
    std::vector<double> dae_dt;
    std::vector<double> dt;
    std::vector<double> tau;
    std::vector<double> flag;
    std::vector<double> problem;
    std::vector<double> solversparam;
    std::vector<double> factor;
    std::vector<double> physicsparam;    
    std::vector<double> externalparam;            
    std::vector<double> vindx;
    std::vector<double> avparam1, avparam2;       
    std::vector<double> stgib;
    std::vector<double> stgdata;
    std::vector<double> stgparam;    
};

struct Mesh {
    std::vector<double> p; // flattened nd × np array, column-major
    std::vector<int> t, tg;// flattened nve × ne array, column-major
    int nd, dim;           // number of spatial dimensions
    int np;                // number of points
    int nve;               // number of vertices per element
    int nvf;               // number of vertices per face
    int nfe;               // number of faces per element
    int ne;                // number of elements    
    int nf;                // number of faces    
    int elemtype;          // elemtype        
    int npe, npf;    
    int nbndexpr, nbcm, nprdexpr, nprdcom;        
    int np_global=0;
    int ne_global=0;
    int nf_global=0;
    
    std::vector<idx_t>  epart_local, elmdist;    

    // For each local element e, elemGlobalID[e] is its global element index
    std::vector<int> elemGlobalID;
    std::vector<int> nodeGlobalID;  // [np]  global node IDs    

    vector<int> f, bf, f2t, t2t, t2f, t2lf, inte, intl, localfaces;    
    std::vector<double> xdg, udg, vdg, wdg, uhat;
    std::vector<int> xdgdims, udgdims, vdgdims, wdgdims, uhatdims, elem2cpu;    
    
    std::vector<int> interfaceConditions;
    std::vector<int> boundaryConditions;
    std::vector<int> curvedBoundaries;
    std::vector<int> periodicBoundaries1;
    std::vector<int> periodicBoundaries2;
    std::vector<int> cartGridPart;
    
    char** boundaryExprs = nullptr;
    char** curvedBoundaryExprs = nullptr;
    char** periodicExprs1 = nullptr;
    char** periodicExprs2 = nullptr;        
};

struct Master 
{    
    vector<double> xpe, gpe, gwe, xpf, gpf, gwf; 
    vector<double> shapeg, shapegt, shapfg, shapfgt, shapent, shapfnt, shapegw, shapfgw, shapen, shapfn;
    vector<double> xp1d, gp1d, gw1d, shap1dg, shap1dgt, shap1dn, shap1dnt, shap1dgw, phielem, phiface;
    vector<int>  telem, tface, perm, permind; 
    int nd, npe, npf, nge, ngf, porder, pgauss, nfe, elemtype, nodetype, nve, nvf, np1d, ng1d, npermind;
};

struct ElementClassification {
    // Local indices of elements on this rank
    std::vector<int> interiorLocal;
    std::vector<int> boundaryLocal;
    std::vector<int> interfaceLocal;

    // Global IDs of those same elements
    std::vector<int> interiorGlobal;
    std::vector<int> boundaryGlobal;
    std::vector<int> interfaceGlobal;

    // Off-rank neighbor elements (unique)
    std::vector<int> neighborElemLocal;
    std::vector<int> neighborElemGlobal; // global ID of neighbor element
    std::vector<int> neighborElemRank;   // owning rank of that neighbor

    // Off-rank outer elements (unique)
    std::vector<int> outerElemLocal;
    std::vector<int> outerElemGlobal; // global ID of outer element
    std::vector<int> outerElemRank;   // owning rank of that outer element    
};

struct DMD {
  std::vector<int> nbsd;                  // neighbors    
  std::vector<int> localelemrecv;         // recv_local_idx
  std::vector<int> localelemsend;         // send_local_idx
  std::vector<std::array<int, 3>> elemrecv; // each row: [sender, recv_local_idx, sender_global_idx]
  std::vector<std::array<int, 3>> elemsend; // each row: [receiver, send_local_idx, recv_global_idx]
  std::vector<int> elempart;              // global element IDs in the partition
  std::vector<int> elempart_local;        // local element IDs in the partition
  std::vector<int> elem2cpu;              // processor ID for each element in the partition
  std::vector<int> elemsendpts;           // number of elements sent to each neighbor
  std::vector<int> elemrecvpts;           // number of elements received from each neighbor
  std::vector<int> elempartpts;           // partition sizes: [interior, interface, exterior]
  std::vector<int> intepartpts;           // optional: [interior, interface1, interface2, exterior]
  std::vector<int> nbinfo;                // neighboring information 
  std::vector<int> bf;
  int numneigh;                           // number of neighbors 
};

#endif