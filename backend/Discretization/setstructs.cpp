/*
    setstructs.cpp

    This file contains functions for initializing and setting up various data structures used in the Exasim backend for high-order finite element/volume methods. The main structures initialized here are commonstruct, appstruct, masterstruct, meshstruct, solstruct, resstruct, and tempstruct. These structures store metadata, mesh information, solution variables, and temporary arrays required for numerical computations.

    Functions:

    - setcommonstruct(commonstruct &common, appstruct &app, masterstruct &master, meshstruct &mesh, string filein, string fileout, Int curvedMesh, Int fileoffset)
        Initializes the commonstruct with parameters from app, master, and mesh structures, as well as input/output file names and mesh curvature information. Sets up MPI-related fields, problem flags, solver parameters, and allocates memory for time integration coefficients and communication buffers.

    - setresstruct(resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, Int backend)
        Allocates and initializes memory for residual arrays (Rq, Ru, Rh) and their sizes in the resstruct. Handles additional allocations if Enzyme AD is enabled.

    - settempstruct(tempstruct &tmp, appstruct &app, masterstruct &master, meshstruct &mesh, Int backend)
        Allocates temporary arrays for intermediate computations, communication buffers for MPI, and sets their sizes in tempstruct. Handles different spatial schemes and backend types.

    - cpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, tempstruct &tmp, commonstruct &common, string filein, string fileout, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank)
        Reads input data, initializes solution, residual, temporary, and common structures for CPU execution. Allocates memory for solution arrays, mesh permutation arrays, and sets up index arrays for communication. Handles mesh curvature and precomputes shape function products.

    - devappstruct(appstruct &dapp, appstruct &app, commonstruct &common)
        Allocates and copies appstruct arrays to device memory for GPU execution.

    - devsolstruct(solstruct &dsol, solstruct &sol, commonstruct &common)
        Allocates and copies solstruct arrays to device memory for GPU execution.

    - devmasterstruct(masterstruct &dmaster, masterstruct &master, commonstruct &common)
        Allocates and copies masterstruct arrays to device memory for GPU execution.

    - devmeshstruct(meshstruct &dmesh, meshstruct &mesh, commonstruct &common)
        Allocates and copies meshstruct arrays to device memory for GPU execution. Handles additional arrays for HDG/EDG schemes.

    - gpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, tempstruct &tmp, commonstruct &common, solstruct &hsol, resstruct &hres, appstruct &happ, masterstruct &hmaster, meshstruct &hmesh, tempstruct &htmp, commonstruct &hcommon)
        Initializes all structures for GPU execution, including device memory allocation and copying from host. Sets up CUDA/HIP handles and allocates solution arrays for source terms and auxiliary variables.

    Notes:
    - The file uses conditional compilation for MPI, CUDA, HIP, and Enzyme AD support.
    - Memory allocation is performed using TemplateMalloc and TemplateCopytoDevice for portability across CPU and GPU backends.
    - The code is designed for parallel execution and supports domain decomposition for distributed memory systems.
*/
#ifndef __SETSTRUCTS
#define __SETSTRUCTS

#include "ismeshcurved.cpp"

void setcommonstruct(commonstruct &common, appstruct &app, masterstruct &master, meshstruct &mesh, 
        string filein, string fileout, Int curvedMesh, Int fileoffset)
{                   
    common.filein = filein;
    common.fileout = fileout;
            
    common.mpiRank = app.comm[0];  // MPI rank      
    common.mpiProcs = app.comm[1]; // number of MPI ranks           
    
    common.exasimpath = trimToSubstringAtLastOccurence(common.exasimpath, "Exasim");     
    if (common.exasimpath == "") {      
      std::filesystem::path cwd = std::filesystem::current_path();
      common.exasimpath = trimToSubstringAtLastOccurence(cwd, "Exasim");            
      if (common.exasimpath == "") common.exasimpath = trimToSubstringAtLastOccurence(common.fileout, "Exasim");       
    }
    if (common.mpiRank==0) std::cout << "exasimpath = "<<common.exasimpath<<std::endl;
    
#ifdef HAVE_ENZYME    
    common.enzyme = 1;
#endif            
    common.read_uh = app.read_uh;
    common.nc = app.ndims[5]; // number of compoments of (u, q)
    common.ncu = app.ndims[6];// number of compoments of (u)        
    common.ncq = app.ndims[7];// number of compoments of (q)
    //common.ncp = app.ndims[8];// number of compoments of (p)    
    common.nco = app.ndims[9];// number of compoments of (o)    
    common.nch = app.ndims[10];// number of compoments of (uhat)
    common.ncx = app.ndims[11];// number of compoments of (xdg)        
    common.nce = app.ndims[12];// number of compoments of (output)        
    common.ncw = app.ndims[13];//number of compoments of (w)
    common.nsca = app.ndims[14];// number of components of scalar fields for visualization
    common.nvec = app.ndims[15];// number of components of vector fields for visualization
    common.nten = app.ndims[16];// number of components of tensor fields for visualization
    common.nsurf = app.ndims[17];// number of components of surface fields for visualization, storage, and QoIs
    common.nvqoi = app.ndims[18];// number of volume quantities of interest (QoIs)    

    common.ncm = 1;//number of components of monitor function    
    if (app.flag[1]==1)
        common.ncs = common.nc;  // wave problem
    else if (app.flag[0]==1)
        common.ncs = common.ncu; // time-dependent problem
    else
        common.ncs = 0; // steady-state problem
        
    common.nd = master.ndims[0];     // spatial dimension    
    common.elemtype = master.ndims[1]; 
    common.nodetype = master.ndims[2]; 
    common.porder = master.ndims[3]; 
    common.pgauss = master.ndims[4]; 
    common.npe = master.ndims[5]; // number of nodes on master element
    common.npf = master.ndims[6]; // number of nodes on master face       
    common.nge = master.ndims[7]; // number of gauss points on master element
    common.ngf = master.ndims[8]; // number of gauss poInts on master face          
    common.np1d = master.ndims[9]; // number of node points on 1D element
    common.ng1d = master.ndims[10]; // number of gauss poInts on 1D element          
    
    common.ne = mesh.ndims[1]; // number of elements in this subdomain 
    common.nf = mesh.ndims[2]; // number of faces in this subdomain 
    common.nv = mesh.ndims[3]; // number of vertices in this subdomain       
    common.nfe = mesh.ndims[4]; // number of faces per element        
    common.nbe = mesh.ndims[5]; // number of blocks for elements 
    common.neb = mesh.ndims[6]; // maximum number of elements per block
    common.nbf = mesh.ndims[7]; // number of blocks for faces   
    common.nfb = mesh.ndims[8]; // maximum number of faces per block     
    
    common.ndof = common.npe*common.ncu*common.ne; // number of degrees of freedom of u
    common.ndofq = common.npe*common.ncq*common.ne; // number of degrees of freedom of q
    common.ndofw = common.npe*common.ncw*common.ne; // number of degrees of freedom of w
    common.ndofuhat = common.npf*common.ncu*common.nf; // number of degrees of freedom of uhat
    common.ndofudg = common.npe*common.nc*common.ne; // number of degrees of freedom of udg    
    common.ndofsdg = common.npe*common.ncs*common.ne; // number of degrees of freedom of sdg    
    common.ndofodg = common.npe*common.nco*common.ne; // number of degrees of freedom of odg               
    common.ndofedg = common.npe*common.nce*common.ne; // number of degrees of freedom of edg               
    common.ndofucg = mesh.nsize[12]-1;
            
    common.ntau = app.nsize[8];
    if (common.ntau>0) common.tau0 = app.tau[0];
		
    common.curvedMesh = curvedMesh;        
    common.fileoffset = fileoffset;

    common.tdep = app.flag[0];      // 0: steady-state; 1: time-dependent;  
    common.wave = app.flag[1];
    common.linearProblem = app.flag[2]; // 0: nonlinear problem;  1: linear problem
    common.debugMode = app.flag[3]; // 1: save data to binary files for debugging
    common.matvecOrder = app.flag[4];        
    common.gmresOrthogMethod = app.flag[5];            
    common.preconditioner = app.flag[6];
    common.precMatrixType = app.flag[7];
    common.ptcMatrixType = app.flag[8];
    common.runmode = app.flag[9];
    common.tdfunc = app.flag[10];
    common.source = app.flag[11]; 
    common.modelnumber = app.flag[12]; 
    common.extFhat = app.flag[13];
    common.extUhat = app.flag[14];
    common.extStab = app.flag[15];
    common.subproblem = app.flag[16];
    
    common.tsteps = app.nsize[4];  // number of time steps          
    common.spatialScheme = app.problem[0];   /* 0: HDG; 1: EDG; 2: IEDG, HEDG */
    common.appname = app.problem[1];   /* 0: Euler; 1: Compressible Navier-Stokes; etc. */    
    common.temporalScheme = app.problem[2];  // 0: DIRK; 1: BDF; 2: ERK
    common.torder = app.problem[3];    /* temporal accuracy order */
    common.tstages = app.problem[4];    /* DIRK stages */    
    common.convStabMethod = app.problem[5];  // Flag for convective stabilization tensor. 0: Constant tau, 1: Lax-Friedrichs; 2: Roe.
    common.diffStabMethod = app.problem[6];  // Flag for diffusive stabilization tensor. 0: No diffusive stabilization.
    common.rotatingFrame = app.problem[7];   // Flag for rotating frame. Options: 0: Velocities are respect to a non-rotating frame. 1: Velocities are respect to a rotating frame.
    common.viscosityModel = app.problem[8];  // Flag for viscosity model. 0: Constant kinematic viscosity; 1: Sutherland's law. 2: Hack for viscosity for compressible HIT case in (Johnsen, 2010)
    common.SGSmodel = app.problem[9];        // Flag for sub-grid scale (SGS) model. Only available for 3D solver.
                                        //  0: No SGS model. 1: Static Smagorinsky/Yoshizawa/Knight model. 
                                        //  2: Static WALE/Yoshizawa/Knight model. 3: Static Vreman/Yoshizawa/Knight model.
                                        //  4: Dynamic Smagorinsky/Yoshizawa/Knight model.        
    common.ALEflag = app.problem[10];   // Flag for Arbitrary Lagrangian-Eulerian (ALE) formulation. 0: No ALE; 1: Translation; 2: Translation + rotation; 3: Translation + rotation + deformation
    common.ncAV = app.problem[11];    // Flag for artificial viscosity. 0: No artificial viscosity; 1: Homogeneous artificial viscosity (C. Nguyen's formulation); 2: Hypersonic homogeneous artificial viscosity (C. Nguyen's formulation)
                                        //                                3: Isotropic artificial viscosity (D. Moro's formulation). 4: Latest version of the model (taking the best of all previous models)
                                        //                                8: Density smoothness sensor (Per's approach)    
    common.linearSolver = app.problem[12];  /* 0: GMRES; 1: CG; etc. */      
    common.nonlinearSolverMaxIter = app.problem[13];                
    common.linearSolverMaxIter = app.problem[14];        
    common.gmresRestart = app.problem[15];    
    common.RBdim = app.problem[16];  
    common.saveSolFreq = app.problem[17];    
    common.saveSolOpt = app.problem[18];    
    common.timestepOffset = app.problem[19];    
    common.stgNmode = app.problem[20];    
    common.saveSolBouFreq = app.problem[21];   
    common.ibs = app.problem[22];   
    common.dae_steps = app.problem[23];  // number of dual time steps      
    common.saveResNorm = app.problem[24];   
    common.AVsmoothingIter = app.problem[25]; //Number of times artificial viscosity is smoothed
    common.frozenAVflag = app.problem[26]; // Flag deciding if artificial viscosity is calculated once per non-linear solve or in every residual evluation
                                           //   0: AV not frozen, evaluated every iteration
                                           //   1: AV frozen, evluated once per solve (default)          
    common.ppdegree = app.problem[27]; // degree of polynomial preconditioner
    common.coupledinterface = app.problem[28]; 
    common.coupledcondition = app.problem[29]; 
    common.coupledboundarycondition = app.problem[30];
    common.AVdistfunction = app.problem[31];
    
    common.RBcurrentdim = 0; // current dimension of the reduced basis space
    common.RBremovedind = 0; // the vector to be removed from the RB space and replaced with new vector
    common.Wcurrentdim = 0; // current dimension of the W space
    
    common.nonlinearSolverTol = app.solversparam[0];    
    common.linearSolverTol = app.solversparam[1];
    common.matvecTol = app.solversparam[2];
    common.PTCparam = app.solversparam[3];
    if (common.tdep==1)
        common.time = app.factor[0];
    common.rampFactor = 1.0;   // Ramp factor for artificial viscosity flux        
    common.dae_alpha = app.factor[1];
    common.dae_beta = app.factor[2];
    common.dae_gamma = app.factor[3];
    common.dae_epsilon = app.factor[4];
    
    common.nstgib = app.nsize[11];
    if (common.nstgib > 0) common.stgib = copyarray(app.stgib,app.nsize[11]); 

    //common.eblks = &mesh.eblks[0]; // element blocks
    //common.fblks = &mesh.fblks[0]; // face blocks        
    //common.dt = &app.dt[0]; // face blocks     
    common.eblks = copyarray(mesh.eblks,mesh.nsize[2]); // element blocks
    common.fblks = copyarray(mesh.fblks,mesh.nsize[3]); // face blocks            
    common.dt = copyarray(app.dt,app.nsize[4]); // timestep sizes       
    common.nvindx = app.nsize[12];
    common.vindx = copyarray(app.vindx,app.nsize[12]); 
    common.dae_dt = copyarray(app.dae_dt,app.nsize[13]); // dual timestep sizes           
    common.szinterfacefluxmap = app.nsize[14];
    common.interfacefluxmap = copyarray(app.interfacefluxmap,app.nsize[14]); 
    common.cartgridpart = copyarray(mesh.cartgridpart,mesh.nsize[25]); 
    common.szcartgridpart = mesh.nsize[25];

    common.boundaryConditions = copyarray(mesh.boundaryConditions, mesh.nsize[27]);
    common.intepartpts = copyarray(mesh.intepartpts, mesh.nsize[28]);
    // if (mesh.nsize[25] > 0) TemplateFree(mesh.cartgridpart, 0);
    // if (mesh.nsize[27] > 0) TemplateFree(mesh.boundaryConditions, 0);
    if (mesh.nsize[28] > 0) {
        if (common.intepartpts[1] > 0) common.isd = 1;
        //TemplateFree(mesh.intepartpts, 0);        
    }
    
    if (common.nvqoi > 0) common.qoivolume = (dstype*) malloc(common.nvqoi*sizeof(dstype));
    if (common.nsurf > 0) common.qoisurface = (dstype*) malloc(common.nsurf*sizeof(dstype));

    common.nf0 = 0;
    for (Int j=0; j<common.nbf; j++) {
        Int f1 = common.fblks[3*j]-1;
        Int f2 = common.fblks[3*j+1];    
        Int ib = common.fblks[3*j+2];    
        if (ib==0) common.nf0 += f2 - f1; // number of interior faces        
    }
    
    Int tstages = common.tstages;  
    if (tstages>0) {
        common.DIRKcoeff_c = (dstype*) malloc(tstages*sizeof(dstype));
        common.DIRKcoeff_d = (dstype*) malloc(tstages*tstages*sizeof(dstype));
        common.DIRKcoeff_t = (dstype*) malloc(tstages*sizeof(dstype));
    }
    if (common.torder>0) {
        common.BDFcoeff_c = (dstype*) malloc((common.torder+1)*sizeof(dstype));
        common.BDFcoeff_t = (dstype*) malloc(sizeof(dstype));
    }    
        
    if (common.mpiProcs>1) { // MPI
        common.ne0 = mesh.elempartpts[0]; // number of interior elements
        common.ne1 = mesh.elempartpts[0]+mesh.elempartpts[1]; // added with number of interface elements
        common.ne2 = mesh.elempartpts[0]+mesh.elempartpts[1]+mesh.elempartpts[2]; // added with number of exterior elements        
        common.ndof1 = common.npe*common.ncu*common.ne1; // number of degrees of freedom
        common.ndofq1 = common.npe*common.ncq*common.ne1; // number of degrees of freedom of q
        common.ndofw1 = common.npe*common.ncw*common.ne1; // number of degrees of freedom of w
        common.ndofudg1 = common.npe*common.nc*common.ne1; // number of degrees of freedom of udg    
        common.ndofsdg1 = common.npe*common.ncs*common.ne1; // number of degrees of freedom of sdg    
        common.ndofodg1 = common.npe*common.nco*common.ne1; // number of degrees of freedom of odg           
        common.ndofedg1 = common.npe*common.nce*common.ne1; // number of degrees of freedom of edg           
        
        common.nbe0 = 0;
        common.nbe1 = 0;
        common.nbe2 = 0;
        common.ncie = 0; // number of coupled interface elements
        for (Int j=0; j<common.nbe; j++) {           
            if (common.eblks[3*j+2] == 0)
                common.nbe0 += 1;
            if (common.eblks[3*j+2] <= 1)
                common.nbe1 += 1;
            if (common.eblks[3*j+2] <= 2)
                common.nbe2 += 1;
            if (common.eblks[3*j+2] == -1)
              common.ncie += common.eblks[3*j+1] - common.eblks[3*j+0] + 1;
        }        
        
        common.nbf0 = mesh.ndims[9]; // number of blocks for faces   
        common.nbf1 = mesh.ndims[10]; // number of blocks for faces   
        
        common.nnbsd = mesh.nsize[4];
        common.nelemsend = mesh.nsize[5];
        common.nelemrecv = mesh.nsize[6];    
        common.nbsd = copyarray(mesh.nbsd,mesh.nsize[4]); 
        common.elemsend = copyarray(mesh.elemsend,mesh.nsize[5]); 
        common.elemrecv = copyarray(mesh.elemrecv,mesh.nsize[6]); 
        common.elemsendpts = copyarray(mesh.elemsendpts,mesh.nsize[7]); 
        common.elemrecvpts = copyarray(mesh.elemrecvpts,mesh.nsize[8]); 
                                
        common.nnbintf = mesh.sznbintf;
        common.nfacesend = mesh.szfacesend;
        common.nfacerecv = mesh.szfacerecv;            
        common.nbintf = copyarray(mesh.nbintf,mesh.sznbintf); 
        common.facesend = copyarray(mesh.facesend,mesh.szfacesend); 
        common.facesendpts = copyarray(mesh.facesendpts,mesh.szfacesendpts); 
        common.facerecv = copyarray(mesh.facerecv,mesh.szfacerecv);         
        common.facerecvpts = copyarray(mesh.facerecvpts,mesh.szfacerecvpts); 
        
#ifdef  HAVE_MPI
        common.requests = (MPI_Request *) malloc( 2*(common.nnbsd + common.nnbintf) * sizeof(MPI_Request) );
        common.statuses = (MPI_Status *) malloc( 2*(common.nnbsd + common.nnbintf) * sizeof(MPI_Status) );     
        if (common.spatialScheme==1) {
          //common.ninterfacefaces = getinterfacefaces(mesh.f2e, common.ne1, common.nf);
          //common.ndofuhatinterface = common.ncu*common.npf*common.ninterfacefaces;
        }
#endif        
    }
    else {
        common.nnbsd = 0;
        common.nelemsend = 0;
        common.nelemrecv = 0;            
        common.nbe0 = common.nbe;
        common.nbe1 = common.nbe;
        common.nbe2 = common.nbe;
        common.ne0 = common.ne;
        common.ne1 = common.ne;
        common.ne2 = common.ne;
        common.ndof1 = common.npe*common.ncu*common.ne; // number of degrees of freedom
        common.ndofq1 = common.npe*common.ncq*common.ne; // number of degrees of freedom of q
        common.ndofw1 = common.npe*common.ncw*common.ne; // number of degrees of freedom of w
        common.ndofudg1 = common.npe*common.nc*common.ne; // number of degrees of freedom of udg    
        common.ndofsdg1 = common.npe*common.ncs*common.ne; // number of degrees of freedom of sdg    
        common.ndofodg1 = common.npe*common.nco*common.ne; // number of degrees of freedom of odg                   
        common.ndofedg1 = common.npe*common.nce*common.ne; // number of degrees of freedom of edg                   
    }            
}

void setresstruct(resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, Int backend)
{
    Int ncu = app.ndims[6];    // number of compoments of (u)
    Int ncq = app.ndims[7];    // number of compoments of (q)
    //Int ncp = app.ndims[8];    // number of compoments of (p)
    Int nd = master.ndims[0];  // spatial dimension    
    Int npe = master.ndims[5]; // number of nodes on master element    
    Int ne = mesh.ndims[1];    // number of elements in this subdomain 
    Int npf = master.ndims[6]; // number of nodes on master face   
    Int nf = mesh.ndims[2];    // number of elements in this subdomain 
    Int nfe = mesh.ndims[4]; // number of faces per element            
   
    TemplateMalloc(&res.Rq, 2*npe*ncu*nd*ne, backend);   
    TemplateMalloc(&res.Ru, npe*ncu*ne, backend);
    TemplateMalloc(&res.Rh, max(npf*max(ncu,ncq)*nf, npf*nfe*ncu*ne), backend);
    res.szRq = 2*npe*ncu*nd*ne;
    res.szRu = npe*ncu*ne;
    res.szRh = max(npf*max(ncu,ncq)*nf, npf*nfe*ncu*ne);
    
#ifdef HAVE_ENZYME
    TemplateMalloc(&res.dRq, 2*npe*ncu*nd*ne, backend);   
    TemplateMalloc(&res.dRu, npe*ncu*ne, backend);
    TemplateMalloc(&res.dRh, npf*max(ncu,ncq)*nf, backend);
#endif
    
    // set pointers    
    res.Rqe = &res.Rq[0];
    res.Rqf = &res.Rq[0];        
    res.Rue = &res.Ru[0];
    res.Ruf = &res.Ru[0];
    
#ifdef HAVE_ENZYME                        
    res.dRqe = &res.dRq[0];
    res.dRqf = &res.dRq[0];        
    res.dRue = &res.dRu[0];
    res.dRuf = &res.dRu[0];
#endif                   
}

void settempstruct(tempstruct &tmp, appstruct &app, masterstruct &master, meshstruct &mesh, Int backend)
{               
    Int nc = app.ndims[5]; // number of compoments of (u, q, p)
    Int ncu = app.ndims[6];// number of compoments of (u)        
    Int ncq = app.ndims[7];    // number of compoments of (q)
    Int nco = app.ndims[9];// number of compoments of (o)    
    Int ncx = app.ndims[11];// number of compoments of (xdg)    
    Int ncw = app.ndims[13];// number of compoments of (u)       
    Int nd = master.ndims[0];     // spatial dimension        
    Int npe = master.ndims[5]; // number of nodes on master element
    Int npf = master.ndims[6]; // number of nodes on master face       
    Int nge = master.ndims[7]; // number of gauss points on master element
    Int ngf = master.ndims[8]; // number of gauss poInts on master face
    Int nfe = mesh.ndims[4]; // number of faces per element                        
    Int neb = mesh.ndims[6]; // maximum number of elements per block
    Int nfb = mesh.ndims[8]; // maximum number of faces per block   
    Int ne = mesh.ndims[1]; // number of elements in this subdomain 
    //Int nf = mesh.ndims[2]; // number of faces in this subdomain 
    Int ndofucg = mesh.nsize[12]-1;; //total DOFs for CG field
    Int spatialScheme = app.problem[0];   /* 0: LDG; 1: HDG */

    Int n0 = max(npe*max(nc+ncw,ncx)*neb, npf*(ncu+2*nc+2*ncw)*nfb);                        
    n0 = 2*max(n0, nco*ne*npe);
    Int n1 = max(ncx+2*nd*nd+1, ncu*nd+ncu+ncw+max(nc,ncu*(nd+1)));
    Int n2 = max(ncx+nd+1+nd*(nd-1), ncu+2*ncu*nd+2*nc+2*ncw);    
    Int n3 = max(nge*n1*neb, ngf*n2*nfb);
    n3 = 2*max(n3, ndofucg);
    
    if (spatialScheme > 0) {
      //Int k1 = npe*ncu*npe*ncu*neb + npe*npf*nfe*ncu*ncu*neb + npe*npf*nfe*ncu*ncu*neb + npf*nfe*npf*nfe*ncu*ncu*neb;
      Int k1 = max(npe*ncu*npe*ncu*neb, npf*nfe*npf*nfe*ncu*ncu*neb);
      Int k2 = npf*nfe*neb*ncu + npf*npf*nfe*neb*ncu*ncu + npf*npf*nfe*neb*ncu*ncq + npf*npf*nfe*neb*ncu*ncu;
      Int k3 = npf*ncu*npf*ncu*nfb;
      Int k4 = nge*(nd+1)*ncu*max(ncu,ncq)*neb;       
      n0 = max(n0, k1);
      n0 = max(n0, k2);      
      n0 = max(n0, k3);  
      n0 = max(n0, k4);  

      int nga = nge*neb;
      k1 = nga*ncu*nd + nga*ncu + nga*nc + nga*ncw + nga*ncu*nd*nc + nga*ncu*nd*ncw + nga*ncu*nc + nga*ncu*ncw + nga*ncw*nc; 
      nga = ngf*neb*nfe; 
      k2 = nga*(ncu + nc + nco + ncw + ncw + ncu*nd + ncu*nd*nc + ncu*nd*ncu + ncu*nd*ncw + ncw*nc); 
      n3 = max(n3, k1);
      n3 = max(n3, k2);
    }
        
    TemplateMalloc(&tmp.tempn, n0+n3, backend); 
    tmp.tempg = &tmp.tempn[n0];
    //TemplateMalloc(&tmp.tempg, n3, backend); 
    tmp.sztempn = n0;
    tmp.sztempg = n3;       
#ifdef HAVE_MPI     
    if (spatialScheme == 0) {
      TemplateMalloc(&tmp.buffsend, max(nc,nco)*npe*mesh.nsize[5], backend);
      TemplateMalloc(&tmp.buffrecv, max(nc,nco)*npe*mesh.nsize[6], backend);
      tmp.szbuffsend = max(nc,nco)*npe*mesh.nsize[5];
      tmp.szbuffrecv = max(nc,nco)*npe*mesh.nsize[6];
    }
    else if (spatialScheme == 1) {
      int m = ncu*npf*nfe; 
      int n = max(nc,nco)*npe; // fix bug here
      TemplateMalloc(&tmp.buffsend, max(m*m + m, n)*mesh.nsize[5], backend);
      TemplateMalloc(&tmp.buffrecv, max(m*m + m, n)*mesh.nsize[6], backend);
      tmp.szbuffsend = max(m*m + m, n)*mesh.nsize[5];
      tmp.szbuffrecv = max(m*m + m, n)*mesh.nsize[6];
      
      TemplateMalloc(&tmp.bufffacesend, app.nsize[14]*npf*mesh.szfacesend, backend);
      TemplateMalloc(&tmp.bufffacerecv, app.nsize[14]*npf*mesh.szfacerecv, backend);
      tmp.szbufffacesend = app.nsize[14]*npf*mesh.szfacesend;
      tmp.szbufffacerecv = app.nsize[14]*npf*mesh.szfacerecv;
    }
#endif
}

void cpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        string filein, string fileout, Int mpiprocs, Int mpirank, Int fileoffset, Int omprank) 
{
     
    if (mpirank==0)
        printf("Reading data from binary files \n");
    readInput(app, master, mesh, sol, filein, mpiprocs, mpirank, fileoffset, omprank);            
    
    if (mpirank==0)
        printf("Finish reading data from binary files \n");
    
//     if (mpiprocs != app.ndims[0]) {
//         if (mpirank==0) {
//             printf("# processors = %d, # subdomains = %d\n", mpiprocs, app.ndims[0]);
//             error("Number of MPI processes is incorrect\n");
//         }
//     }
         
    // offset facecon
//     for (Int i=0; i<mesh.nsize[1]; i++)
//         mesh.facecon[i] = mesh.facecon[i] - 1;
                
//     if (mpiprocs > 1) {
//         // offset nbsd
//         for (Int i=0; i<mesh.nsize[4]; i++)
//             mesh.nbsd[i] = mesh.nbsd[i] - 1;
// 
//         // offset elemsend
//         for (Int i=0; i<mesh.nsize[5]; i++)
//             mesh.elemsend[i] = mesh.elemsend[i] - 1;
// 
//         // offset elemrecv
//         for (Int i=0; i<mesh.nsize[6]; i++)
//             mesh.elemrecv[i] = mesh.elemrecv[i] - 1;
//     }
            
    if (mpirank==0) printf("Set res struct... \n");    
    setresstruct(res, app, master, mesh, 0);
    
    if (mpirank==0) printf("Set temp struct... \n");    
    settempstruct(tmp, app, master, mesh, 0);    
        
    if (mpirank==0) printf("Run IsMeshCurved... \n");    
    int curvedmesh = IsMeshCurved(sol, app, master, mesh, tmp);       
    if (mpirank==0) printf("IsMeshCurved = %d \n",curvedmesh);        
    
    if (mpirank==0) printf("Set common struct... \n");
    setcommonstruct(common, app, master, mesh, filein, fileout, curvedmesh, fileoffset);            
    common.cublasHandle = 0;
    common.eventHandle = 0; 
    if (mpirank==0) printf("Finish setting common struct... \n");
                        
    if (common.spatialScheme > 0) {      
      int nd = common.nd;
      int npe = common.npe;
      int nge = common.nge;
      master.shapegwdotshapeg = (dstype*) malloc (sizeof (dstype)*npe*npe*nge*(nd+1));      
      master.szshapegwdotshapeg = npe*npe*nge*(nd+1);
      for (int d=0; d<nd+1; d++) {
        for (int g=0; g<nge; g++) {
          for (int i=0; i<npe; i++) {
            for (int j=0; j<npe; j++) {            
              master.shapegwdotshapeg[j+i*npe+g*npe*npe+d*npe*npe*nge] = master.shapegw[j+g*npe+d*npe*nge]*master.shapegt[g+i*nge];              
            }            
          }
        }
      }      
      int npf = common.npf;
      int ngf = common.ngf;
      master.shapfgwdotshapfg = (dstype*) malloc (sizeof (dstype)*npf*npf*ngf*nd);
      master.szshapfgwdotshapfg = npf*npf*ngf*nd;
      for (int d=0; d<nd; d++) {
        for (int g=0; g<ngf; g++) {
          for (int i=0; i<npf; i++) {
            for (int j=0; j<npf; j++) {            
              master.shapfgwdotshapfg[j+i*npf+g*npf*npf+d*npf*npf*ngf] = master.shapfgw[j+g*npf+d*npf*ngf]*master.shapfgt[g+i*ngf]; // fix bug here             
            }            
          }
        }
      }              
      // if (common.nelemsend>0) {
      //     mesh.interfacefaces = (Int*) malloc (sizeof (Int)*common.nelemsend);          
      //     int n = getsubdomaininterfaces(mesh.interfacefaces, mesh.f2e, common.ne1, common.nf);
      //     if (n != common.nelemsend) error("Number of interfaces mismatch");
      // }
      
      Int nf = common.nf;
      Int ne = common.ne;
      Int nfe = common.nfe;
      mesh.e2f = (Int*) malloc (sizeof (Int)*nfe*ne);
      mesh.f2f = (Int*) malloc (sizeof (Int)*2*(nfe-1)*nf);
      mesh.f2l = (Int*) malloc (sizeof (Int)*2*(nfe-1)*nf);
      mke2f(mesh.e2f, mesh.f2e, nf, nfe, ne);
      mkf2f(mesh.f2f, mesh.f2l, mesh.f2e, mesh.e2f, nf, nfe, ne);       
      
//       for (int i = 0; i < common.nelemsend; i++) {
//         int esend = common.elemsend[i];
//         int erecv = common.elemrecv[i];
//         for (int j = 0; j < nfe; j++) 
//           for (int k = 0; k < nfe; k++) {
//             if (mesh.e2f[j + nfe*esend] == mesh.e2f[k + nfe*erecv]) {
//               int f = mesh.e2f[j + nfe*esend];
//               printf("%d %d %d %d %d %d %d %d %d %d %d ", common.mpiRank, i, esend, j, erecv, k, f, mesh.f2e[0 + 4*f], mesh.f2e[1 + 4*f], mesh.f2e[2 + 4*f], mesh.f2e[3 + 4*f]);
//               for (int l = 0; l< 2*(nfe-1); l++)
//                 printf("%d ", mesh.f2f[l + 2*(nfe-1)*f]);
//               printf("\n");
//             }
//         }
//       }
      
//       int n = getsubdomaininterfaces(mesh.f2e, common.ne1, common.nf);
//       int *interfacefaces = (Int*) malloc (sizeof (Int)*n);          
//       getsubdomaininterfaces(interfacefaces, mesh.f2e, common.ne1, common.nf);
//       print2iarray(interfacefaces, 1, n);
//       print2iarray(mesh.f2e, 4, n);
//       print2iarray(mesh.f2f, 2*(nfe-1), n);
//       print2iarray(mesh.f2l, 2*(nfe-1), n);
//       CPUFREE(interfacefaces);      
      
//       print2iarray(mesh.f2e, 4, nf);
//       print2iarray(mesh.e2f, nfe, ne);
//       print2iarray(mesh.f2f, 2*(nfe-1), nf);
//       print2iarray(mesh.f2l, 2*(nfe-1), nf);
//      CPUFREE(e2f);      
    }

    if (common.ncs>0) {
        // initialize source term
        Int N = common.npe*common.ncs*common.ne;
        sol.sdg = (dstype*) malloc (sizeof (dstype)*N);        
        for (Int i=0; i<N; i++)
            sol.sdg[i] = 0.0;               
        sol.sdgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.ncs*common.ne);     
        sol.szsdg = N;
        sol.szsdgg = common.nge*common.ncs*common.ne;
    }            
    
    if (common.ncw>0) {        
        Int N = common.npe*common.ncw*common.ne;
        sol.wsrc = (dstype*) malloc (sizeof (dstype)*N);
        sol.szwsrc = N;
        for (Int i=0; i<N; i++)
            sol.wsrc[i] = 0.0;               
        if (common.dae_steps>0) {
            sol.wdual = (dstype*) malloc (sizeof (dstype)*N);
            sol.szwdual = N;
        }
    }
    
    if (common.compudgavg>0) {
        sol.udgavg = (dstype*) malloc (sizeof (dstype)*(common.npe*common.nc*common.ne1+1));  
        sol.szudgavg = common.npe*common.nc*common.ne1+1;
    }
    
    // allocate memory for uh
    if (!app.read_uh) {
    sol.uh = (dstype*) malloc (sizeof (dstype)*common.npf*common.ncu*common.nf);
    }
    sol.szuh = common.npf*common.ncu*common.nf;
    
    #ifdef HAVE_ENZYME
        sol.duh = (dstype*) malloc (sizeof (dstype)*common.npf*common.ncu*common.nf);
    #endif
    //sol.uhg = (dstype*) malloc (sizeof (dstype)*common.ngf*common.ncu*common.nf);
    //sol.udgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nc*common.ne);
    if (common.nco>0) {
        sol.odgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
        sol.og1 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
        sol.og2 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
        sol.szodgg = common.nge*common.nco*common.ne;
        sol.szog1 = common.ngf*common.nco*common.nf;
        sol.szog2 = common.ngf*common.nco*common.nf;
        
        #ifdef HAVE_ENZYME
            sol.dodgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
            ArraySetValue(sol.dodgg, zero, common.nge*common.nco*common.ne, 0);
            sol.dog1 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
            ArraySetValue(sol.dog1, zero, common.ngf*common.nco*common.nf, 0);
            sol.dog2 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
            ArraySetValue(sol.dog2, zero, common.ngf*common.nco*common.nf, 0);
        #endif
    }
            
    if (mpirank==0) printf("Precompute index arrays... \n");    
    //mesh.index = (Int*) malloc (sizeof (Int)*(1024));            
    
    // allocate memory 
    mesh.findxdg1 = (Int*) malloc (sizeof (Int)*common.npf*common.ncx*common.nf);
    mesh.findxdgp = (Int*) malloc (sizeof (Int)*(common.nbf+1));
    mesh.findudg1 = (Int*) malloc (sizeof (Int)*common.npf*common.nc*common.nf);
    mesh.findudg2 = (Int*) malloc (sizeof (Int)*common.npf*common.nc*common.nf);
    mesh.findudgp = (Int*) malloc (sizeof (Int)*(common.nbf+1));
    mesh.eindudg1 = (Int*) malloc (sizeof (Int)*common.npe*common.nc*common.ne);
    mesh.eindudgp = (Int*) malloc (sizeof (Int)*(common.nbe+1));
    
    mesh.szfindxdg1 = common.npf*common.ncx*common.nf;
    mesh.szfindxdgp = (common.nbf+1);
    mesh.szfindudg1 = common.npf * common.nc * common.nf;
    mesh.szfindudg2 = common.npf * common.nc * common.nf;
    mesh.szfindudgp = common.nbf + 1;
    mesh.szeindudg1 = common.npe * common.nc * common.ne;
    mesh.szeindudgp = common.nbe + 1;
            
    faceperm1(mesh.findxdg1, mesh.findxdgp, mesh.facecon, mesh.fblks,
             common.npf, common.ncx, common.npe, common.ncx, common.nbf);
    faceperm(mesh.findudg1, mesh.findudg2, mesh.findudgp, mesh.facecon, mesh.fblks,
             common.npf, common.nc, common.npe, common.nc, common.nbf);
    elemperm(mesh.eindudg1, mesh.eindudgp, mesh.eblks, 
             common.npe, common.nc, common.nc, common.nbe);   
    
    if (common.nelemsend>0) {
        Int bsz = common.npe*common.ncu;
        Int nudg = common.npe*common.nc;
        mesh.elemsendind = (Int*) malloc (sizeof (Int)*bsz*common.nelemsend);                        
        mesh.szelemsendind = bsz*common.nelemsend;
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++) 
                mesh.elemsendind[bsz*i+j] = nudg*common.elemsend[i] + j;            
      
        bsz = common.npe*common.ncAV;
        nudg = common.npe*common.nco;
        mesh.elemsendodg = (Int*) malloc (sizeof (Int)*bsz*common.nelemsend);
        mesh.szelemsendodg = bsz*common.nelemsend;
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemsendodg[bsz*i+j] = nudg*common.elemsend[i] + j;
       
        bsz = common.npe*common.nc;
        nudg = common.npe*common.nc;
        mesh.elemsendudg = (Int*) malloc (sizeof (Int)*bsz*common.nelemsend);
        mesh.szelemsendudg = bsz*common.nelemsend;
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemsendudg[bsz*i+j] = nudg*common.elemsend[i] + j;

    }
    
    if (common.nelemrecv>0) {
        Int bsz = common.npe*common.ncu;
        Int nudg = common.npe*common.nc;
        mesh.elemrecvind = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);                        
        mesh.szelemrecvind = bsz*common.nelemrecv;
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++) 
                mesh.elemrecvind[bsz*i+j] = nudg*common.elemrecv[i] + j;            
    
        bsz = common.npe*common.ncAV;
        nudg = common.npe*common.nco;
        mesh.elemrecvodg = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);
        mesh.szelemrecvodg = bsz*common.nelemrecv;
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemrecvodg[bsz*i+j] = nudg*common.elemrecv[i] + j;  
    
        bsz = common.npe*common.nc;
        nudg = common.npe*common.nc;
        mesh.elemrecvudg = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);
        mesh.szelemrecvudg = bsz*common.nelemrecv;
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemrecvudg[bsz*i+j] = nudg*common.elemrecv[i] + j;
    }            
    
    if (mpirank==0) {
        //print2darray(app.physicsparam,1,app.nsize[6]); 
        printf("finish cpuInit... \n");
    }
}

#ifdef HAVE_GPU

void devappstruct(appstruct &dapp, appstruct &app, commonstruct &common)
{        
    TemplateMalloc(&dapp.nsize, app.lsize[0], common.backend);
    TemplateMalloc(&dapp.ndims, app.nsize[0], common.backend);
    TemplateMalloc(&dapp.flag, app.nsize[1], common.backend);    
    TemplateMalloc(&dapp.problem, app.nsize[2], common.backend);    
    TemplateMalloc(&dapp.uinf, app.nsize[3], common.backend);    
    TemplateMalloc(&dapp.dt, app.nsize[4], common.backend);    
    TemplateMalloc(&dapp.factor, app.nsize[5], common.backend);    
    TemplateMalloc(&dapp.physicsparam, app.nsize[6], common.backend);    
    TemplateMalloc(&dapp.solversparam, app.nsize[7], common.backend);  
    TemplateMalloc(&dapp.tau, app.nsize[8], common.backend);  
    TemplateMalloc(&dapp.stgdata, app.nsize[9], common.backend);  
    TemplateMalloc(&dapp.stgparam, app.nsize[10], common.backend);  
    TemplateMalloc(&dapp.stgib, app.nsize[11], common.backend);  
    TemplateMalloc(&dapp.vindx, app.nsize[12], common.backend);  
    TemplateMalloc(&dapp.dae_dt, app.nsize[13], common.backend);  
    TemplateMalloc(&dapp.interfacefluxmap, app.nsize[14], common.backend);  
    TemplateMalloc(&dapp.avparam, app.nsize[15], common.backend);  
    
    TemplateCopytoDevice( dapp.nsize, app.nsize, app.lsize[0], common.backend );      
    TemplateCopytoDevice( dapp.ndims, app.ndims, app.nsize[0], common.backend );      
    TemplateCopytoDevice( dapp.flag, app.flag, app.nsize[1], common.backend );      
    TemplateCopytoDevice( dapp.problem, app.problem, app.nsize[2], common.backend );      
    TemplateCopytoDevice( dapp.uinf, app.uinf, app.nsize[3], common.backend );      
    TemplateCopytoDevice( dapp.dt, app.dt, app.nsize[4], common.backend );      
    TemplateCopytoDevice( dapp.factor, app.factor, app.nsize[5], common.backend );      
    TemplateCopytoDevice( dapp.physicsparam, app.physicsparam, app.nsize[6], common.backend );      
    TemplateCopytoDevice( dapp.solversparam, app.solversparam, app.nsize[7], common.backend );          
    TemplateCopytoDevice( dapp.tau, app.tau, app.nsize[8], common.backend );          
    TemplateCopytoDevice( dapp.stgdata, app.stgdata, app.nsize[9], common.backend );          
    TemplateCopytoDevice( dapp.stgparam, app.stgparam, app.nsize[10], common.backend );          
    TemplateCopytoDevice( dapp.stgib, app.stgib, app.nsize[11], common.backend );   
    TemplateCopytoDevice( dapp.vindx, app.vindx, app.nsize[12], common.backend );   
    TemplateCopytoDevice( dapp.dae_dt, app.dae_dt, app.nsize[13], common.backend );   
    TemplateCopytoDevice( dapp.interfacefluxmap, app.interfacefluxmap, app.nsize[14], common.backend );   
    TemplateCopytoDevice( dapp.avparam, app.avparam, app.nsize[15], common.backend );   
    
    dapp.szflag = app.nsize[1];
    dapp.szproblem = app.nsize[2];
    dapp.szuinf = app.nsize[3];
    dapp.szdt = app.nsize[4];
    dapp.szfactor = app.nsize[5];
    dapp.szphysicsparam = app.nsize[6];
    dapp.szsolversparam = app.nsize[7];
    dapp.sztau = app.nsize[8];
    dapp.szstgdata = app.nsize[9];
    dapp.szstgparam = app.nsize[10];
    dapp.szstgib = app.nsize[11];
    dapp.szvindx = app.nsize[12];
    dapp.szdae_dt = app.nsize[13];
    dapp.szinterfacefluxmap = app.nsize[14];
    dapp.szavparam = app.nsize[15];
    
    Int ncu, ncq, ncw;
    ncu = app.ndims[6];// number of compoments of (u)
    ncq = app.ndims[7];// number of compoments of (q)
    ncw = app.ndims[13];// number of compoments of (w)    
    if (ncu>0) {
        TemplateMalloc(&dapp.fc_u, ncu, common.backend); 
        TemplateCopytoDevice( dapp.fc_u, app.fc_u, ncu, common.backend );          
        TemplateMalloc(&dapp.dtcoef_u, ncu, common.backend); 
        TemplateCopytoDevice( dapp.dtcoef_u, app.dtcoef_u, ncu, common.backend );          
    }        
    if (ncq>0) {
        TemplateMalloc(&dapp.fc_q, ncq, common.backend); 
        TemplateCopytoDevice( dapp.fc_q, app.fc_q, ncq, common.backend );           
        TemplateMalloc(&dapp.dtcoef_q, ncq, common.backend); 
        TemplateCopytoDevice( dapp.dtcoef_q, app.dtcoef_q, ncq, common.backend );          
    }        
    if (ncw>0) {
        TemplateMalloc(&dapp.fc_w, ncw, common.backend); 
        TemplateCopytoDevice( dapp.fc_w, app.fc_w, ncw, common.backend );          
        TemplateMalloc(&dapp.dtcoef_w, ncw, common.backend); 
        TemplateCopytoDevice( dapp.dtcoef_w, app.dtcoef_w, ncw, common.backend );          
    }                    
}

void devsolstruct(solstruct &dsol, solstruct &sol, commonstruct &common)
{    
    TemplateMalloc(&dsol.nsize, sol.lsize[0], common.backend);
    TemplateMalloc(&dsol.ndims, sol.nsize[0], common.backend);
    TemplateMalloc(&dsol.xdg, sol.nsize[1], common.backend);    
    TemplateMalloc(&dsol.udg, sol.nsize[2], common.backend);    
    //TemplateMalloc(&dsol.uh, sol.nsize[3], common.backend);    
    TemplateMalloc(&dsol.odg, sol.nsize[3], common.backend);          
    TemplateMalloc(&dsol.wdg, sol.nsize[4], common.backend);

    #ifdef HAVE_ENZYME
    TemplateMalloc(&dsol.dudg, sol.nsize[2], common.backend);    
    //TemplateMalloc(&dsol.uh, sol.nsize[3], common.backend);    
    TemplateMalloc(&dsol.dodg, sol.nsize[3], common.backend);          
    TemplateMalloc(&dsol.dwdg, sol.nsize[4], common.backend);
    #endif          
    
    TemplateCopytoDevice( dsol.nsize, sol.nsize, sol.lsize[0], common.backend );          
    TemplateCopytoDevice( dsol.ndims, sol.ndims, sol.nsize[0], common.backend );      
    TemplateCopytoDevice( dsol.xdg, sol.xdg, sol.nsize[1], common.backend );      
    TemplateCopytoDevice( dsol.udg, sol.udg, sol.nsize[2], common.backend );      
    //TemplateCopytoDevice( dsol.uh, sol.uh, sol.nsize[3], common.backend );      
    TemplateCopytoDevice( dsol.odg, sol.odg, sol.nsize[3], common.backend );      
    TemplateCopytoDevice( dsol.wdg, sol.wdg, sol.nsize[4], common.backend );   

    #ifdef HAVE_ENZYME
    TemplateCopytoDevice( dsol.dudg, sol.dudg, sol.nsize[2], common.backend );   
    //CHECK( cudaMemcpy( dsol.uh, sol.uh, sol.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    TemplateCopytoDevice( dsol.dodg, sol.dodg, sol.nsize[3], common.backend );   
    TemplateCopytoDevice( dsol.wdg, sol.dwdg, sol.nsize[4], common.backend );   
    #endif
}

void devmasterstruct(masterstruct &dmaster, masterstruct &master, commonstruct &common)
{    
    TemplateMalloc(&dmaster.nsize, master.lsize[0], common.backend);
    TemplateMalloc(&dmaster.ndims, master.nsize[0], common.backend);
    TemplateMalloc(&dmaster.shapegt, master.nsize[1], common.backend);    
    TemplateMalloc(&dmaster.shapegw, master.nsize[2], common.backend);    
    TemplateMalloc(&dmaster.shapfgt, master.nsize[3], common.backend);    
    TemplateMalloc(&dmaster.shapfgw, master.nsize[4], common.backend);    
    TemplateMalloc(&dmaster.shapent, master.nsize[5], common.backend);    
    TemplateMalloc(&dmaster.shapen, master.nsize[6], common.backend);    
    TemplateMalloc(&dmaster.shapfnt, master.nsize[7], common.backend);    
    TemplateMalloc(&dmaster.shapfn, master.nsize[8], common.backend);        
    TemplateMalloc(&dmaster.xpe, master.nsize[9], common.backend);    
    TemplateMalloc(&dmaster.gpe, master.nsize[10], common.backend);    
    TemplateMalloc(&dmaster.gwe, master.nsize[11], common.backend);    
    TemplateMalloc(&dmaster.xpf, master.nsize[12], common.backend);    
    TemplateMalloc(&dmaster.gpf, master.nsize[13], common.backend);    
    TemplateMalloc(&dmaster.gwf, master.nsize[14], common.backend);    
    TemplateMalloc(&dmaster.shap1dgt, master.nsize[15], common.backend);    
    TemplateMalloc(&dmaster.shap1dgw, master.nsize[16], common.backend);    
    TemplateMalloc(&dmaster.shap1dnt, master.nsize[17], common.backend);    
    TemplateMalloc(&dmaster.shap1dnl, master.nsize[18], common.backend);    
    TemplateMalloc(&dmaster.xp1d, master.nsize[19], common.backend);    
    TemplateMalloc(&dmaster.gp1d, master.nsize[20], common.backend);    
    TemplateMalloc(&dmaster.gw1d, master.nsize[21], common.backend);    
    
    TemplateCopytoDevice( dmaster.nsize, master.nsize, master.lsize[0], common.backend);    
    TemplateCopytoDevice( dmaster.ndims, master.ndims, master.nsize[0], common.backend );      
    TemplateCopytoDevice( dmaster.shapegt, master.shapegt, master.nsize[1], common.backend );      
    TemplateCopytoDevice( dmaster.shapegw, master.shapegw, master.nsize[2], common.backend );      
    TemplateCopytoDevice( dmaster.shapfgt, master.shapfgt, master.nsize[3], common.backend );      
    TemplateCopytoDevice( dmaster.shapfgw, master.shapfgw, master.nsize[4], common.backend );      
    TemplateCopytoDevice( dmaster.shapent, master.shapent, master.nsize[5], common.backend );      
    TemplateCopytoDevice( dmaster.shapen, master.shapen, master.nsize[6], common.backend );      
    TemplateCopytoDevice( dmaster.shapfnt, master.shapfnt, master.nsize[7], common.backend );      
    TemplateCopytoDevice( dmaster.shapfn, master.shapfn, master.nsize[8], common.backend );      
    TemplateCopytoDevice( dmaster.xpe, master.xpe, master.nsize[9], common.backend );      
    TemplateCopytoDevice( dmaster.gpe, master.gpe, master.nsize[10], common.backend );      
    TemplateCopytoDevice( dmaster.gwe, master.gwe, master.nsize[11], common.backend );      
    TemplateCopytoDevice( dmaster.xpf, master.xpf, master.nsize[12], common.backend );      
    TemplateCopytoDevice( dmaster.gpf, master.gpf, master.nsize[13], common.backend );      
    TemplateCopytoDevice( dmaster.gwf, master.gwf, master.nsize[14], common.backend );          
    TemplateCopytoDevice( dmaster.shap1dgt, master.shap1dgt, master.nsize[15], common.backend );          
    TemplateCopytoDevice( dmaster.shap1dgw, master.shap1dgw, master.nsize[16], common.backend );          
    TemplateCopytoDevice( dmaster.shap1dnt, master.shap1dnt, master.nsize[17], common.backend );          
    TemplateCopytoDevice( dmaster.shap1dnl, master.shap1dnl, master.nsize[18], common.backend );          
    TemplateCopytoDevice( dmaster.xp1d, master.xp1d, master.nsize[19], common.backend );          
    TemplateCopytoDevice( dmaster.gp1d, master.gp1d, master.nsize[20], common.backend );          
    TemplateCopytoDevice( dmaster.gw1d, master.gw1d, master.nsize[21], common.backend );          
}

void devmeshstruct(meshstruct &dmesh, meshstruct &mesh, commonstruct &common)
{
    TemplateMalloc(&dmesh.nsize, mesh.lsize[0], common.backend);
    TemplateMalloc(&dmesh.ndims, mesh.nsize[0], common.backend);
    TemplateMalloc(&dmesh.facecon, mesh.nsize[1], common.backend);
    TemplateMalloc(&dmesh.eblks, mesh.nsize[2], common.backend);
    TemplateMalloc(&dmesh.fblks, mesh.nsize[3], common.backend);
    TemplateMalloc(&dmesh.nbsd, mesh.nsize[4], common.backend);
    TemplateMalloc(&dmesh.elemsend, mesh.nsize[5], common.backend);
    TemplateMalloc(&dmesh.elemrecv, mesh.nsize[6], common.backend);
    TemplateMalloc(&dmesh.elemsendpts, mesh.nsize[7], common.backend);
    TemplateMalloc(&dmesh.elemrecvpts, mesh.nsize[8], common.backend);
    TemplateMalloc(&dmesh.elempart, mesh.nsize[9], common.backend);
    TemplateMalloc(&dmesh.elempartpts, mesh.nsize[10], common.backend);
    TemplateMalloc(&dmesh.cgelcon, mesh.nsize[11], common.backend);
    TemplateMalloc(&dmesh.rowent2elem, mesh.nsize[12], common.backend);
    TemplateMalloc(&dmesh.cgent2dgent, mesh.nsize[13], common.backend);
    TemplateMalloc(&dmesh.colent2elem, mesh.nsize[14], common.backend);
    TemplateMalloc(&dmesh.rowe2f1, mesh.nsize[15], common.backend);
    TemplateMalloc(&dmesh.cole2f1, mesh.nsize[16], common.backend);
    TemplateMalloc(&dmesh.ent2ind1, mesh.nsize[17], common.backend);
    TemplateMalloc(&dmesh.rowe2f2, mesh.nsize[18], common.backend);
    TemplateMalloc(&dmesh.cole2f2, mesh.nsize[19], common.backend);
    TemplateMalloc(&dmesh.ent2ind2, mesh.nsize[20], common.backend);

    TemplateCopytoDevice( dmesh.nsize, mesh.nsize, mesh.lsize[0], common.backend);
    TemplateCopytoDevice(dmesh.ndims, mesh.ndims, mesh.nsize[0], common.backend);
    TemplateCopytoDevice(dmesh.facecon, mesh.facecon, mesh.nsize[1], common.backend);
    TemplateCopytoDevice(dmesh.eblks, mesh.eblks, mesh.nsize[2], common.backend);
    TemplateCopytoDevice(dmesh.fblks, mesh.fblks, mesh.nsize[3], common.backend);
    TemplateCopytoDevice(dmesh.nbsd, mesh.nbsd, mesh.nsize[4], common.backend);
    TemplateCopytoDevice(dmesh.elemsend, mesh.elemsend, mesh.nsize[5], common.backend);
    TemplateCopytoDevice(dmesh.elemrecv, mesh.elemrecv, mesh.nsize[6], common.backend);
    TemplateCopytoDevice(dmesh.elemsendpts, mesh.elemsendpts, mesh.nsize[7], common.backend);
    TemplateCopytoDevice(dmesh.elemrecvpts, mesh.elemrecvpts, mesh.nsize[8], common.backend);
    TemplateCopytoDevice(dmesh.elempart, mesh.elempart, mesh.nsize[9], common.backend);
    TemplateCopytoDevice(dmesh.elempartpts, mesh.elempartpts, mesh.nsize[10], common.backend);
    TemplateCopytoDevice(dmesh.cgelcon, mesh.cgelcon, mesh.nsize[11], common.backend);
    TemplateCopytoDevice(dmesh.rowent2elem, mesh.rowent2elem, mesh.nsize[12], common.backend);
    TemplateCopytoDevice(dmesh.cgent2dgent, mesh.cgent2dgent, mesh.nsize[13], common.backend);
    TemplateCopytoDevice(dmesh.colent2elem, mesh.colent2elem, mesh.nsize[14], common.backend);
    TemplateCopytoDevice(dmesh.rowe2f1, mesh.rowe2f1, mesh.nsize[15], common.backend);
    TemplateCopytoDevice(dmesh.cole2f1, mesh.cole2f1, mesh.nsize[16], common.backend);
    TemplateCopytoDevice(dmesh.ent2ind1, mesh.ent2ind1, mesh.nsize[17], common.backend);
    TemplateCopytoDevice(dmesh.rowe2f2, mesh.rowe2f2, mesh.nsize[18], common.backend);
    TemplateCopytoDevice(dmesh.cole2f2, mesh.cole2f2, mesh.nsize[19], common.backend);
    TemplateCopytoDevice(dmesh.ent2ind2, mesh.ent2ind2, mesh.nsize[20], common.backend);
    
    if (common.spatialScheme > 0) {      
        TemplateMalloc(&dmesh.f2e, mesh.nsize[21], common.backend);
        TemplateMalloc(&dmesh.elemcon, mesh.nsize[22], common.backend);
        TemplateMalloc(&dmesh.perm, mesh.nsize[23], common.backend);
        TemplateCopytoDevice(dmesh.f2e, mesh.f2e, mesh.nsize[21], common.backend);
        TemplateCopytoDevice(dmesh.elemcon, mesh.elemcon, mesh.nsize[22], common.backend);
        TemplateCopytoDevice(dmesh.perm, mesh.perm, mesh.nsize[23], common.backend);
        
        if (mesh.szfaceperm>0) {
          TemplateMalloc(&dmesh.faceperm, mesh.szfaceperm, common.backend);
          TemplateCopytoDevice(dmesh.faceperm, mesh.faceperm, mesh.szfaceperm, common.backend);
        }          
        if (mesh.sznbintf>0) {
          TemplateMalloc(&dmesh.nbintf, mesh.sznbintf, common.backend);
          TemplateCopytoDevice(dmesh.nbintf, mesh.nbintf, mesh.sznbintf, common.backend);
        }          
        if (mesh.szfacesend>0) {
          TemplateMalloc(&dmesh.facesend, mesh.szfacesend, common.backend);
          TemplateCopytoDevice(dmesh.facesend, mesh.facesend, mesh.szfacesend, common.backend);
        }          
        if (mesh.szfacesendpts>0) {
          TemplateMalloc(&dmesh.facesendpts, mesh.szfacesendpts, common.backend);
          TemplateCopytoDevice(dmesh.facesendpts, mesh.facesendpts, mesh.szfacesendpts, common.backend);
        }          
        if (mesh.szfacerecv>0) {
          TemplateMalloc(&dmesh.facerecv, mesh.szfacerecv, common.backend);
          TemplateCopytoDevice(dmesh.facerecv, mesh.facerecv, mesh.szfacerecv, common.backend);
        }          
        if (mesh.szfacerecvpts>0) {
          TemplateMalloc(&dmesh.facerecvpts, mesh.szfacerecvpts, common.backend);
          TemplateCopytoDevice(dmesh.facerecvpts, mesh.facerecvpts, mesh.szfacerecvpts, common.backend);
        }          
    }

    //cudaTemplateMalloc(&dmesh.index, 1024);
    
    Int nbe = mesh.ndims[5];
    Int nbf = mesh.ndims[7];
    TemplateMalloc(&dmesh.findxdg1, mesh.findxdgp[nbf], common.backend);
    TemplateMalloc(&dmesh.findxdgp, nbf+1, common.backend);
    TemplateMalloc(&dmesh.findudg1, mesh.findudgp[nbf], common.backend);
    TemplateMalloc(&dmesh.findudg2, mesh.findudgp[nbf], common.backend);
    TemplateMalloc(&dmesh.findudgp, nbf+1, common.backend);
    TemplateMalloc(&dmesh.eindudg1, mesh.eindudgp[nbe], common.backend);
    TemplateMalloc(&dmesh.eindudgp, nbe+1, common.backend);
            
    TemplateCopytoDevice(dmesh.findxdg1, mesh.findxdg1, mesh.findxdgp[nbf], common.backend);
    TemplateCopytoDevice(dmesh.findxdgp, mesh.findxdgp, nbf+1, common.backend);
    TemplateCopytoDevice(dmesh.findudg1, mesh.findudg1, mesh.findudgp[nbf], common.backend);
    TemplateCopytoDevice(dmesh.findudg2, mesh.findudg2, mesh.findudgp[nbf], common.backend);
    TemplateCopytoDevice(dmesh.findudgp, mesh.findudgp, nbf+1, common.backend);
    TemplateCopytoDevice(dmesh.eindudg1, mesh.eindudg1, mesh.eindudgp[nbe], common.backend);
    TemplateCopytoDevice(dmesh.eindudgp, mesh.eindudgp, nbe+1, common.backend);
        
    if (common.nelemsend>0) {
        Int bsz = common.npe*common.ncu*common.nelemsend;
        TemplateMalloc(&dmesh.elemsendind, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemsendind, mesh.elemsendind, bsz, common.backend);
    
        bsz = common.npe*common.nc*common.nelemsend;
        TemplateMalloc(&dmesh.elemsendudg, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemsendudg, mesh.elemsendudg, bsz, common.backend);
    
        bsz = common.npe*common.ncAV*common.nelemsend;
        TemplateMalloc(&dmesh.elemsendodg, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemsendodg, mesh.elemsendodg, bsz, common.backend);
    }
    
    if (common.nelemrecv>0) {
        Int bsz = common.npe*common.ncu*common.nelemrecv;
        TemplateMalloc(&dmesh.elemrecvind, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemrecvind, mesh.elemrecvind, bsz, common.backend);
    
        bsz = common.npe*common.nc*common.nelemrecv;
        TemplateMalloc(&dmesh.elemrecvudg, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemrecvudg, mesh.elemrecvudg, bsz, common.backend);
     
        bsz = common.npe*common.ncAV*common.nelemrecv;
        TemplateMalloc(&dmesh.elemrecvodg, bsz, common.backend);
        TemplateCopytoDevice(dmesh.elemrecvodg, mesh.elemrecvodg, bsz, common.backend);
    }
}

void gpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
       meshstruct &mesh, tempstruct &tmp, commonstruct &common, solstruct &hsol, resstruct &hres, 
       appstruct &happ, masterstruct &hmaster, meshstruct &hmesh, tempstruct &htmp, commonstruct &hcommon) 
{    
    devappstruct(app, happ, hcommon);
    devmasterstruct(master, hmaster, hcommon);    
    devmeshstruct(mesh, hmesh, hcommon);
    devsolstruct(sol, hsol, hcommon);    
    setresstruct(res, happ, hmaster, hmesh, hcommon.backend);    
    settempstruct(tmp, happ, hmaster, hmesh, hcommon.backend);            
                
    // set common struct
    setcommonstruct(common, happ, hmaster, hmesh, 
            hcommon.filein, hcommon.fileout, hcommon.curvedMesh, hcommon.fileoffset);        
    
    int mpirank = hcommon.mpiRank;
    if (mpirank==0) printf("Finish setting common struct... \n");

    if (common.spatialScheme > 0) {
//       if (common.nelemsend > 0) {
//         TemplateMalloc(&mesh.interfacefaces, common.nelemsend, common.backend);       
//         TemplateCopytoDevice(mesh.interfacefaces, hmesh.interfacefaces, common.nelemsend, common.backend);          
//       }

      int M = common.npe*common.npe*common.nge*(common.nd+1);
      TemplateMalloc(&master.shapegwdotshapeg, M, common.backend);       
      TemplateCopytoDevice(master.shapegwdotshapeg, hmaster.shapegwdotshapeg, M, common.backend);

      M = common.npf*common.npf*common.ngf*(common.nd);
      TemplateMalloc(&master.shapfgwdotshapfg, M, common.backend);       
      TemplateCopytoDevice(master.shapfgwdotshapfg, hmaster.shapfgwdotshapfg, M, common.backend);        
      
      M = 2*(common.nfe-1)*common.nf;
      TemplateMalloc(&mesh.f2f, M, common.backend);       
      TemplateMalloc(&mesh.f2l, M, common.backend);       
      TemplateCopytoDevice(mesh.f2f, hmesh.f2f, M, common.backend);
      TemplateCopytoDevice(mesh.f2l, hmesh.f2l, M, common.backend);
      TemplateMalloc(&mesh.e2f, common.nfe*common.ne, common.backend);
      TemplateCopytoDevice(mesh.e2f, hmesh.e2f, common.nfe*common.ne, common.backend);
    }

#ifdef HAVE_CUDA    
    // create cuda event handle
    CHECK(cudaEventCreate(&common.eventHandle));
    
    // create cublas handle
    CHECK_CUBLAS(cublasCreate(&common.cublasHandle));
    CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_HOST));                     //     CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_DEVICE));    
#endif        

#ifdef HAVE_HIP    
    // create cuda event handle
    CHECK(hipEventCreate(&common.eventHandle));
    
    // create cublas handle
    CHECK_HIPBLAS(hipblasCreate(&common.cublasHandle));
    CHECK_HIPBLAS(hipblasSetPointerMode(common.cublasHandle, HIPBLAS_POINTER_MODE_HOST));                     //     CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_DEVICE));    
#endif        
    
    if (common.ncs>0) {
        // initialize source term
        Int N = common.npe*common.ncs*common.ne;
        TemplateMalloc(&sol.sdg, N, common.backend);       
        TemplateCopytoDevice(sol.sdg, hsol.sdg, N, common.backend);
        TemplateMalloc(&sol.sdgg, common.nge*common.ncs*common.ne, common.backend);              
    }          
    
    if (common.ncw>0) {        
        Int N = common.npe*common.ncw*common.ne;
        TemplateMalloc(&sol.wsrc, N, common.backend);       
        TemplateCopytoDevice(sol.wsrc, hsol.wsrc, N, common.backend);
        if (common.dae_steps>0)
            TemplateMalloc(&sol.wdual, N, common.backend);                 
    }          
    
    if (common.compudgavg>0) {
        TemplateMalloc(&sol.udgavg, common.npe*common.nc*common.ne1+1, common.backend);
    }
    
    TemplateMalloc(&sol.uh, common.npf*common.ncu*common.nf, common.backend);    
    if (common.read_uh) {
        TemplateCopytoDevice(sol.uh, hsol.uh, common.npf*common.ncu*common.nf, common.backend);
    }

    #ifdef HAVE_ENZYME
        TemplateMalloc(&sol.duh, common.npf*common.ncu*common.nf, common.backend);
    #endif
    //cudaTemplateMalloc(&sol.uhg, common.ngf*common.ncu*common.nf);   
    //cudaTemplateMalloc(&sol.udgg, common.nge*common.nc*common.ne);   
    if (common.nco>0) {
        TemplateMalloc(&sol.odgg, common.nge*common.nco*common.ne, common.backend);
        TemplateMalloc(&sol.og1, common.ngf*common.nco*common.nf, common.backend);    
        TemplateMalloc(&sol.og2, common.ngf*common.nco*common.nf, common.backend);    
        #ifdef HAVE_ENZYME
        TemplateMalloc(&sol.dodgg, common.nge*common.nco*common.ne, common.backend);
        TemplateCopytoDevice(sol.dodgg, hsol.dodgg, common.nge*common.nco*common.ne, common.backend);

        TemplateMalloc(&sol.dog1, common.ngf*common.nco*common.nf, common.backend);    
        TemplateCopytoDevice(sol.dog1, hsol.dog1, common.ngf*common.nco*common.nf, common.backend);

        TemplateMalloc(&sol.dog2, common.ngf*common.nco*common.nf, common.backend); 
        TemplateCopytoDevice(sol.dog2, hsol.dog2, common.ngf*common.nco*common.nf, common.backend);
        #endif
    }
    
    if (common.mpiRank==0) 
        printf("finish gpuInit... \n");    
}
#endif

#endif
