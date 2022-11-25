#ifndef __SETSTRUCTS
#define __SETSTRUCTS

#include "ismeshcurved.cpp"

void setcommonstruct(commonstruct &common, appstruct &app, masterstruct &master, meshstruct &mesh, 
        string filein, string fileout, Int curvedMesh)
{               
    common.filein = filein;
    common.fileout = fileout;
            
    common.mpiRank = app.comm[0];  // MPI rank      
    common.mpiProcs = app.comm[1]; // number of MPI ranks           
    
#ifdef HAVE_ENZYME    
    common.enzyme = 1;
#endif            
    
    common.nc = app.ndims[5]; // number of compoments of (u, q)
    common.ncu = app.ndims[6];// number of compoments of (u)        
    common.ncq = app.ndims[7];// number of compoments of (q)
    //common.ncp = app.ndims[8];// number of compoments of (p)    
    common.nco = app.ndims[9];// number of compoments of (o)    
    common.nch = app.ndims[10];// number of compoments of (uhat)
    common.ncx = app.ndims[11];// number of compoments of (xdg)        
    common.nce = app.ndims[12];// number of compoments of (output)        
    common.ncw = app.ndims[13];//number of compoments of (w)    
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
    
    common.ndof = common.npe*common.ncu*common.ne; // number of degrees of freedom
    common.ndofq = common.npe*common.ncq*common.ne; // number of degrees of freedom of q
    common.ndofw = common.npe*common.ncw*common.ne; // number of degrees of freedom of w
    common.ndofudg = common.npe*common.nc*common.ne; // number of degrees of freedom of udg    
    common.ndofsdg = common.npe*common.ncs*common.ne; // number of degrees of freedom of sdg    
    common.ndofodg = common.npe*common.nco*common.ne; // number of degrees of freedom of odg               
    common.ndofedg = common.npe*common.nce*common.ne; // number of degrees of freedom of edg               
    common.ndofucg = mesh.nsize[12]-1;
            
    common.ntau = app.nsize[8];
    if (common.ntau>0)
	common.tau0 = app.tau[0];
		
    common.curvedMesh = curvedMesh;        

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
    common.stgib = copyarray(app.stgib,app.nsize[11]); 

    //common.eblks = &mesh.eblks[0]; // element blocks
    //common.fblks = &mesh.fblks[0]; // face blocks        
    //common.dt = &app.dt[0]; // face blocks     
    common.eblks = copyarray(mesh.eblks,mesh.nsize[2]); // element blocks
    common.fblks = copyarray(mesh.fblks,mesh.nsize[3]); // face blocks            
    common.dt = copyarray(app.dt,app.nsize[4]); // timestep sizes       
    common.nvindx = app.nsize[12];
    common.vindx = copyarray(app.vindx,app.nsize[12]); 
    common.dae_dt = copyarray(app.dae_dt,app.nsize[13]); // dual timestep sizes   
    
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
        for (Int j=0; j<common.nbe; j++) {           
            if (common.eblks[3*j+2] <= 0)
                common.nbe0 += 1;
            if (common.eblks[3*j+2] <= 1)
                common.nbe1 += 1;
            if (common.eblks[3*j+2] <= 2)
                common.nbe2 += 1;
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
        
#ifdef  HAVE_MPI
        common.requests = (MPI_Request *) malloc( 2*common.nnbsd * sizeof(MPI_Request) );
        common.statuses = (MPI_Status *) malloc( 2*common.nnbsd * sizeof(MPI_Status) );     
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

void setresstruct(resstruct &res, appstruct &app, masterstruct &master, meshstruct &mesh, Int cpuMemory)
{
    Int ncu = app.ndims[6];    // number of compoments of (u)
    Int ncq = app.ndims[7];    // number of compoments of (q)
    //Int ncp = app.ndims[8];    // number of compoments of (p)
    Int nd = master.ndims[0];  // spatial dimension    
    Int npe = master.ndims[5]; // number of nodes on master element    
    Int ne = mesh.ndims[1];    // number of elements in this subdomain 
    Int npf = master.ndims[6]; // number of nodes on master face   
    Int nf = mesh.ndims[2];    // number of elements in this subdomain 
    
    if (cpuMemory==0) { // GPU memory
#ifdef HAVE_CUDA    
        cudaTemplateMalloc(&res.Rq, 2*npe*ncu*nd*ne);    
        cudaTemplateMalloc(&res.Ru, npe*ncu*ne); 
        cudaTemplateMalloc(&res.Rh, npf*max(ncu,ncq)*nf); 
#ifdef HAVE_ENZYME                    
        cudaTemplateMalloc(&res.dRq, 2*npe*ncu*nd*ne);    
        cudaTemplateMalloc(&res.dRu, npe*ncu*ne); 
        cudaTemplateMalloc(&res.dRh, npf*max(ncu,ncq)*nf);         
#endif                    
#endif            
    }
    else {// CPU memory
        res.Rq = (dstype *) malloc(2*npe*ncu*nd*ne*sizeof(dstype));
        res.Ru = (dstype *) malloc((npe*ncu*ne)*sizeof(dstype));
        res.Rh = (dstype *) malloc((npf*max(ncu,ncq)*nf)*sizeof(dstype));
#ifdef HAVE_ENZYME                    
        res.dRq = (dstype *) malloc(2*npe*ncu*nd*ne*sizeof(dstype));
        res.dRu = (dstype *) malloc((npe*ncu*ne)*sizeof(dstype));
        res.dRh = (dstype *) malloc((npf*max(ncu,ncq)*nf)*sizeof(dstype));
#endif                            
    }
    
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

void settempstruct(tempstruct &tmp, appstruct &app, masterstruct &master, meshstruct &mesh, Int cpuMemory)
{               
    Int nc = app.ndims[5]; // number of compoments of (u, q, p)
    Int ncu = app.ndims[6];// number of compoments of (u)        
    Int nco = app.ndims[9];// number of compoments of (o)    
    Int ncx = app.ndims[11];// number of compoments of (xdg)    
    Int ncw = app.ndims[13];// number of compoments of (u)       
    Int nd = master.ndims[0];     // spatial dimension        
    Int npe = master.ndims[5]; // number of nodes on master element
    Int npf = master.ndims[6]; // number of nodes on master face       
    Int nge = master.ndims[7]; // number of gauss points on master element
    Int ngf = master.ndims[8]; // number of gauss poInts on master face              
    Int neb = mesh.ndims[6]; // maximum number of elements per block
    Int nfb = mesh.ndims[8]; // maximum number of faces per block   

    Int ne = mesh.ndims[1]; // number of elements in this subdomain 

    Int ndofucg = mesh.nsize[12]-1;; //total DOFs for CG field

    Int n0 = max(npe*max(nc+ncw,ncx)*neb, npf*(ncu+2*nc+2*ncw)*nfb);                        
    n0 = max(n0, nco*ne*npe);
    Int n1 = max(ncx+2*nd*nd+1, ncu*nd+ncu+ncw+max(nc,ncu*(nd+1)));
    Int n2 = max(ncx+nd+1+nd*(nd-1), ncu+2*ncu*nd+2*nc+2*ncw);    
    Int n3 = max(nge*n1*neb, ngf*n2*nfb);
    n3 = max(n3, ndofucg);
    
    #ifdef HAVE_ENZYME
        n0 = 2*n0;
        n3 = 2*n3;
    #endif

    if (cpuMemory==0) { // GPU memory
#ifdef HAVE_CUDA    
        cudaTemplateMalloc(&tmp.tempn, n0);  
        cudaTemplateMalloc(&tmp.tempg, n3);       
#ifdef  HAVE_MPI    
        cudaTemplateMalloc(&tmp.buffsend, max(nc,nco)*npe*mesh.nsize[5]);  
        cudaTemplateMalloc(&tmp.buffrecv, max(nc,nco)*npe*mesh.nsize[6]);       
#endif                      
#endif            
    }
    else { // CPU memory    
        tmp.tempn = (dstype *) malloc(n0*sizeof(dstype));
        tmp.tempg = (dstype *) malloc(n3*sizeof(dstype));
#ifdef  HAVE_MPI            
        tmp.buffsend = (dstype *) malloc(max(nc,nco)*npe*mesh.nsize[5]*sizeof(dstype));
        tmp.buffrecv = (dstype *) malloc(max(nc,nco)*npe*mesh.nsize[6]*sizeof(dstype));
#endif                    
    }        
}

void cpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
        meshstruct &mesh, tempstruct &tmp, commonstruct &common,
        string filein, string fileout, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank) 
{
    if (mpirank==0)
        printf("Reading data from binary files \n");
    readInput(app, master, mesh, sol, filein, mpiprocs, mpirank, ompthreads, omprank);            
    
    if (mpiprocs != app.ndims[0]) {
        if (mpirank==0) 
            error("Number of MPI processes is incorrect\n");
    }
         
    // offset facecon
    for (Int i=0; i<mesh.nsize[1]; i++)
        mesh.facecon[i] = mesh.facecon[i] - 1;
                
    if (mpiprocs > 1) {
        // offset nbsd
        for (Int i=0; i<mesh.nsize[4]; i++)
            mesh.nbsd[i] = mesh.nbsd[i] - 1;

        // offset elemsend
        for (Int i=0; i<mesh.nsize[5]; i++)
            mesh.elemsend[i] = mesh.elemsend[i] - 1;

        // offset elemrecv
        for (Int i=0; i<mesh.nsize[6]; i++)
            mesh.elemrecv[i] = mesh.elemrecv[i] - 1;
    }
            
    setresstruct(res, app, master, mesh, 1);
    settempstruct(tmp, app, master, mesh, 1);    
                
    int curvedmesh = IsMeshCurved(sol, app, master, mesh, tmp);   
        
    if (mpirank==0) 
        printf("Set common struct... \n");
    setcommonstruct(common, app, master, mesh, filein, fileout, curvedmesh);        
    common.cpuMemory = 1;
    common.cublasHandle = 0;
    common.eventHandle = 0;                             
        
    if (common.ncs>0) {
        // initialize source term
        Int N = common.npe*common.ncs*common.ne;
        sol.sdg = (dstype*) malloc (sizeof (dstype)*N);
        for (Int i=0; i<N; i++)
            sol.sdg[i] = 0.0;               
        sol.sdgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.ncs*common.ne);     
    }            
    
    if (common.ncw>0) {        
        Int N = common.npe*common.ncw*common.ne;
        sol.wsrc = (dstype*) malloc (sizeof (dstype)*N);
        for (Int i=0; i<N; i++)
            sol.wsrc[i] = 0.0;               
        if (common.dae_steps>0)
            sol.wdual = (dstype*) malloc (sizeof (dstype)*N);
    }
    
    if (common.compudgavg>0)
        sol.udgavg = (dstype*) malloc (sizeof (dstype)*(common.npe*common.nc*common.ne1+1));  
    
    // allocate memory for uh
    sol.uh = (dstype*) malloc (sizeof (dstype)*common.npf*common.ncu*common.nf);
    #ifdef HAVE_ENZYME
        sol.duh = (dstype*) malloc (sizeof (dstype)*common.npf*common.ncu*common.nf);
    #endif
    //sol.uhg = (dstype*) malloc (sizeof (dstype)*common.ngf*common.ncu*common.nf);
    //sol.udgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nc*common.ne);
    if (common.nco>0) {
        sol.odgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
        sol.og1 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
        sol.og2 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
        #ifdef HAVE_ENZYME
            sol.dodgg = (dstype*) malloc (sizeof (dstype)*common.nge*common.nco*common.ne);
            ArraySetValue(sol.dodgg, zero, common.nge*common.nco*common.ne, 0);
            sol.dog1 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
            ArraySetValue(sol.dog1, zero, common.ngf*common.nco*common.nf, 0);
            sol.dog2 = (dstype*) malloc (sizeof (dstype)*common.ngf*common.nco*common.nf);
            ArraySetValue(sol.dog2, zero, common.ngf*common.nco*common.nf, 0);
        #endif
    }
    
    if (mpirank==0) 
        printf("Precompute index arrays... \n");    
    mesh.index = (Int*) malloc (sizeof (Int)*(1024));            
    
    // allocate memory 
    mesh.findxdg1 = (Int*) malloc (sizeof (Int)*common.npf*common.ncx*common.nf);
    mesh.findxdgp = (Int*) malloc (sizeof (Int)*(common.nbf+1));
    mesh.findudg1 = (Int*) malloc (sizeof (Int)*common.npf*common.nc*common.nf);
    mesh.findudg2 = (Int*) malloc (sizeof (Int)*common.npf*common.nc*common.nf);
    mesh.findudgp = (Int*) malloc (sizeof (Int)*(common.nbf+1));
    mesh.eindudg1 = (Int*) malloc (sizeof (Int)*common.npe*common.nc*common.ne);
    mesh.eindudgp = (Int*) malloc (sizeof (Int)*(common.nbe+1));
    
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
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++) 
                mesh.elemsendind[bsz*i+j] = nudg*common.elemsend[i] + j;            
      
        bsz = common.npe*common.ncAV;
        nudg = common.npe*common.nco;
        mesh.elemsendodg = (Int*) malloc (sizeof (Int)*bsz*common.nelemsend);
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemsendodg[bsz*i+j] = nudg*common.elemsend[i] + j;
       
        bsz = common.npe*common.nc;
        nudg = common.npe*common.nc;
        mesh.elemsendudg = (Int*) malloc (sizeof (Int)*bsz*common.nelemsend);
        for (Int i=0; i<common.nelemsend; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemsendudg[bsz*i+j] = nudg*common.elemsend[i] + j;

    }
    
    if (common.nelemrecv>0) {
        Int bsz = common.npe*common.ncu;
        Int nudg = common.npe*common.nc;
        mesh.elemrecvind = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);                        
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++) 
                mesh.elemrecvind[bsz*i+j] = nudg*common.elemrecv[i] + j;            
    
        bsz = common.npe*common.ncAV;
        nudg = common.npe*common.nco;
        mesh.elemrecvodg = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemrecvodg[bsz*i+j] = nudg*common.elemrecv[i] + j;  
    
        bsz = common.npe*common.nc;
        nudg = common.npe*common.nc;
        mesh.elemrecvudg = (Int*) malloc (sizeof (Int)*bsz*common.nelemrecv);
        for (Int i=0; i<common.nelemrecv; i++)
            for (Int j=0; j<bsz; j++)
                mesh.elemrecvudg[bsz*i+j] = nudg*common.elemrecv[i] + j;
    }            
    
    if (mpirank==0) {
        //print2darray(app.physicsparam,1,app.nsize[6]); 
        printf("finish cpuInit... \n");
    }
}

#ifdef HAVE_CUDA

void devappstruct(appstruct &dapp, appstruct &app)
{        
    cudaTemplateMalloc(&dapp.nsize, app.lsize[0]);
    cudaTemplateMalloc(&dapp.ndims, app.nsize[0]);
    cudaTemplateMalloc(&dapp.flag, app.nsize[1]);    
    cudaTemplateMalloc(&dapp.problem, app.nsize[2]);    
    cudaTemplateMalloc(&dapp.uinf, app.nsize[3]);    
    cudaTemplateMalloc(&dapp.dt, app.nsize[4]);    
    cudaTemplateMalloc(&dapp.factor, app.nsize[5]);    
    cudaTemplateMalloc(&dapp.physicsparam, app.nsize[6]);    
    cudaTemplateMalloc(&dapp.solversparam, app.nsize[7]);  
    cudaTemplateMalloc(&dapp.tau, app.nsize[8]);  
    cudaTemplateMalloc(&dapp.stgdata, app.nsize[9]);  
    cudaTemplateMalloc(&dapp.stgparam, app.nsize[10]);  
    
    CHECK( cudaMemcpy( dapp.nsize, app.nsize, app.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.ndims, app.ndims, app.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.flag, app.flag, app.nsize[1]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.problem, app.problem, app.nsize[2]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.uinf, app.uinf, app.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.dt, app.dt, app.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.factor, app.factor, app.nsize[5]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.physicsparam, app.physicsparam, app.nsize[6]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dapp.solversparam, app.solversparam, app.nsize[7]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dapp.tau, app.tau, app.nsize[8]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dapp.stgdata, app.stgdata, app.nsize[9]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dapp.stgparam, app.stgparam, app.nsize[10]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    
    Int ncu, ncq, ncp;
    ncu = app.ndims[6];// number of compoments of (u)
    ncq = app.ndims[7];// number of compoments of (q)
    ncp = app.ndims[8];// number of compoments of (p)    
    if (ncu>0) {
        cudaTemplateMalloc(&dapp.fc_u, ncu); 
        CHECK( cudaMemcpy( dapp.fc_u, app.fc_u, ncu*sizeof(dstype), cudaMemcpyHostToDevice ) );  
        cudaTemplateMalloc(&dapp.dtcoef_u, ncu); 
        CHECK( cudaMemcpy( dapp.dtcoef_u, app.dtcoef_u, ncu*sizeof(dstype), cudaMemcpyHostToDevice ) );  
    }        
    if (ncq>0) {
        cudaTemplateMalloc(&dapp.fc_q, ncq); 
        CHECK( cudaMemcpy( dapp.fc_q, app.fc_q, ncq*sizeof(dstype), cudaMemcpyHostToDevice ) );       
        cudaTemplateMalloc(&dapp.dtcoef_q, ncq); 
        CHECK( cudaMemcpy( dapp.dtcoef_q, app.dtcoef_q, ncq*sizeof(dstype), cudaMemcpyHostToDevice ) );  
    }        
    if (ncp>0) {
        cudaTemplateMalloc(&dapp.fc_p, ncp); 
        CHECK( cudaMemcpy( dapp.fc_p, app.fc_p, ncp*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        cudaTemplateMalloc(&dapp.dtcoef_p, ncp); 
        CHECK( cudaMemcpy( dapp.dtcoef_p, app.dtcoef_p, ncp*sizeof(dstype), cudaMemcpyHostToDevice ) );              
    }                    
}


void devsolstruct(solstruct &dsol, solstruct &sol)
{    
    cudaTemplateMalloc(&dsol.nsize, sol.lsize[0]);
    cudaTemplateMalloc(&dsol.ndims, sol.nsize[0]);
    cudaTemplateMalloc(&dsol.xdg, sol.nsize[1]);    
    cudaTemplateMalloc(&dsol.udg, sol.nsize[2]);    
    //cudaTemplateMalloc(&dsol.uh, sol.nsize[3]);    
    cudaTemplateMalloc(&dsol.odg, sol.nsize[3]);          
    cudaTemplateMalloc(&dsol.wdg, sol.nsize[4]);

    #ifdef HAVE_ENZYME
    cudaTemplateMalloc(&dsol.dudg, sol.nsize[2]);    
    //cudaTemplateMalloc(&dsol.uh, sol.nsize[3]);    
    cudaTemplateMalloc(&dsol.dodg, sol.nsize[3]);          
    cudaTemplateMalloc(&dsol.dwdg, sol.nsize[4]);
    #endif          
    
    CHECK( cudaMemcpy( dsol.nsize, sol.nsize, sol.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.ndims, sol.ndims, sol.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.xdg, sol.xdg, sol.nsize[1]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.udg, sol.udg, sol.nsize[2]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    //CHECK( cudaMemcpy( dsol.uh, sol.uh, sol.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.odg, sol.odg, sol.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.wdg, sol.wdg, sol.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );   

    #ifdef HAVE_ENZYME
    CHECK( cudaMemcpy( dsol.dudg, sol.dudg, sol.nsize[2]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    //CHECK( cudaMemcpy( dsol.uh, sol.uh, sol.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.dodg, sol.dodg, sol.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dsol.wdg, sol.dwdg, sol.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );  
    #endif
}

void devmasterstruct(masterstruct &dmaster, masterstruct &master)
{    
    cudaTemplateMalloc(&dmaster.nsize, master.lsize[0]);
    cudaTemplateMalloc(&dmaster.ndims, master.nsize[0]);
    cudaTemplateMalloc(&dmaster.shapegt, master.nsize[1]);    
    cudaTemplateMalloc(&dmaster.shapegw, master.nsize[2]);    
    cudaTemplateMalloc(&dmaster.shapfgt, master.nsize[3]);    
    cudaTemplateMalloc(&dmaster.shapfgw, master.nsize[4]);    
    cudaTemplateMalloc(&dmaster.shapent, master.nsize[5]);    
    cudaTemplateMalloc(&dmaster.shapen, master.nsize[6]);    
    cudaTemplateMalloc(&dmaster.shapfnt, master.nsize[7]);    
    cudaTemplateMalloc(&dmaster.shapfn, master.nsize[8]);        
    cudaTemplateMalloc(&dmaster.xpe, master.nsize[9]);    
    cudaTemplateMalloc(&dmaster.gpe, master.nsize[10]);    
    cudaTemplateMalloc(&dmaster.gwe, master.nsize[11]);    
    cudaTemplateMalloc(&dmaster.xpf, master.nsize[12]);    
    cudaTemplateMalloc(&dmaster.gpf, master.nsize[13]);    
    cudaTemplateMalloc(&dmaster.gwf, master.nsize[14]);    
    cudaTemplateMalloc(&dmaster.shap1dgt, master.nsize[15]);    
    cudaTemplateMalloc(&dmaster.shap1dgw, master.nsize[16]);    
    cudaTemplateMalloc(&dmaster.shap1dnt, master.nsize[17]);    
    cudaTemplateMalloc(&dmaster.shap1dnl, master.nsize[18]);    
    cudaTemplateMalloc(&dmaster.xp1d, master.nsize[19]);    
    cudaTemplateMalloc(&dmaster.gp1d, master.nsize[20]);    
    cudaTemplateMalloc(&dmaster.gw1d, master.nsize[21]);    
    
    CHECK( cudaMemcpy( dmaster.nsize, master.nsize, master.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.ndims, master.ndims, master.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapegt, master.shapegt, master.nsize[1]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapegw, master.shapegw, master.nsize[2]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapfgt, master.shapfgt, master.nsize[3]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapfgw, master.shapfgw, master.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapent, master.shapent, master.nsize[5]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapen, master.shapen, master.nsize[6]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapfnt, master.shapfnt, master.nsize[7]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.shapfn, master.shapfn, master.nsize[8]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.xpe, master.xpe, master.nsize[9]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.gpe, master.gpe, master.nsize[10]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.gwe, master.gwe, master.nsize[11]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.xpf, master.xpf, master.nsize[12]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.gpf, master.gpf, master.nsize[13]*sizeof(dstype), cudaMemcpyHostToDevice ) );      
    CHECK( cudaMemcpy( dmaster.gwf, master.gwf, master.nsize[14]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.shap1dgt, master.shap1dgt, master.nsize[15]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.shap1dgw, master.shap1dgw, master.nsize[16]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.shap1dnt, master.shap1dnt, master.nsize[17]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.shap1dnl, master.shap1dnl, master.nsize[18]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.xp1d, master.xp1d, master.nsize[19]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.gp1d, master.gp1d, master.nsize[20]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
    CHECK( cudaMemcpy( dmaster.gw1d, master.gw1d, master.nsize[21]*sizeof(dstype), cudaMemcpyHostToDevice ) );          
}

void devmeshstruct(meshstruct &dmesh, meshstruct &mesh, commonstruct &common)
{
    cudaTemplateMalloc(&dmesh.nsize, mesh.lsize[0]);
    cudaTemplateMalloc(&dmesh.ndims, mesh.nsize[0]);
    cudaTemplateMalloc(&dmesh.facecon, mesh.nsize[1]);
    cudaTemplateMalloc(&dmesh.eblks, mesh.nsize[2]);
    cudaTemplateMalloc(&dmesh.fblks, mesh.nsize[3]);
    cudaTemplateMalloc(&dmesh.nbsd, mesh.nsize[4]);
    cudaTemplateMalloc(&dmesh.elemsend, mesh.nsize[5]);
    cudaTemplateMalloc(&dmesh.elemrecv, mesh.nsize[6]);
    cudaTemplateMalloc(&dmesh.elemsendpts, mesh.nsize[7]);
    cudaTemplateMalloc(&dmesh.elemrecvpts, mesh.nsize[8]);
    cudaTemplateMalloc(&dmesh.elempart, mesh.nsize[9]);
    cudaTemplateMalloc(&dmesh.elempartpts, mesh.nsize[10]);
    cudaTemplateMalloc(&dmesh.cgelcon, mesh.nsize[11]);
    cudaTemplateMalloc(&dmesh.rowent2elem, mesh.nsize[12]);
    cudaTemplateMalloc(&dmesh.cgent2dgent, mesh.nsize[13]);
    cudaTemplateMalloc(&dmesh.colent2elem, mesh.nsize[14]);
    cudaTemplateMalloc(&dmesh.rowe2f1, mesh.nsize[15]);
    cudaTemplateMalloc(&dmesh.cole2f1, mesh.nsize[16]);
    cudaTemplateMalloc(&dmesh.ent2ind1, mesh.nsize[17]);
    cudaTemplateMalloc(&dmesh.rowe2f2, mesh.nsize[18]);
    cudaTemplateMalloc(&dmesh.cole2f2, mesh.nsize[19]);
    cudaTemplateMalloc(&dmesh.ent2ind2, mesh.nsize[20]);
    
    CHECK( cudaMemcpy( dmesh.nsize, mesh.nsize, mesh.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.ndims, mesh.ndims, mesh.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.facecon, mesh.facecon, mesh.nsize[1]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.eblks, mesh.eblks, mesh.nsize[2]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.fblks, mesh.fblks, mesh.nsize[3]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.nbsd, mesh.nbsd, mesh.nsize[4]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.elemsend, mesh.elemsend, mesh.nsize[5]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.elemrecv, mesh.elemrecv, mesh.nsize[6]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.elemsendpts, mesh.elemsendpts, mesh.nsize[7]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.elemrecvpts, mesh.elemrecvpts, mesh.nsize[8]*sizeof(Int), cudaMemcpyHostToDevice ) );    
    CHECK( cudaMemcpy( dmesh.elempart, mesh.elempart, mesh.nsize[9]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.elempartpts, mesh.elempartpts, mesh.nsize[10]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.cgelcon, mesh.cgelcon, mesh.nsize[11]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.rowent2elem, mesh.rowent2elem, mesh.nsize[12]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.cgent2dgent, mesh.cgent2dgent, mesh.nsize[13]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.colent2elem, mesh.colent2elem, mesh.nsize[14]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.rowe2f1, mesh.rowe2f1, mesh.nsize[15]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.cole2f1, mesh.cole2f1, mesh.nsize[16]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.ent2ind1, mesh.ent2ind1, mesh.nsize[17]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.rowe2f2, mesh.rowe2f2, mesh.nsize[18]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.cole2f2, mesh.cole2f2, mesh.nsize[19]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.ent2ind2, mesh.ent2ind2, mesh.nsize[20]*sizeof(Int), cudaMemcpyHostToDevice ) );
    
    cudaTemplateMalloc(&dmesh.index, 1024);
    
    Int nbe = mesh.ndims[5];
    Int nbf = mesh.ndims[7];
    cudaTemplateMalloc(&dmesh.findxdg1, mesh.findxdgp[nbf]);
    cudaTemplateMalloc(&dmesh.findxdgp, nbf+1);
    cudaTemplateMalloc(&dmesh.findudg1, mesh.findudgp[nbf]);
    cudaTemplateMalloc(&dmesh.findudg2, mesh.findudgp[nbf]);
    cudaTemplateMalloc(&dmesh.findudgp, nbf+1);
    cudaTemplateMalloc(&dmesh.eindudg1, mesh.eindudgp[nbe]);
    cudaTemplateMalloc(&dmesh.eindudgp, nbe+1);
            
    CHECK( cudaMemcpy( dmesh.findxdg1, mesh.findxdg1, mesh.findxdgp[nbf]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.findxdgp, mesh.findxdgp, (nbf+1)*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.findudg1, mesh.findudg1, mesh.findudgp[nbf]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.findudg2, mesh.findudg2, mesh.findudgp[nbf]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.findudgp, mesh.findudgp, (nbf+1)*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.eindudg1, mesh.eindudg1, mesh.eindudgp[nbe]*sizeof(Int), cudaMemcpyHostToDevice ) );
    CHECK( cudaMemcpy( dmesh.eindudgp, mesh.eindudgp, (nbe+1)*sizeof(Int), cudaMemcpyHostToDevice ) );
        
    if (common.nelemsend>0) {
        Int bsz = common.npe*common.ncu*common.nelemsend;
        cudaTemplateMalloc(&dmesh.elemsendind, bsz);
        CHECK( cudaMemcpy( dmesh.elemsendind , mesh.elemsendind , bsz*sizeof(Int), cudaMemcpyHostToDevice ) );        
    
        bsz = common.npe*common.nc*common.nelemsend;
        cudaTemplateMalloc(&dmesh.elemsendudg, bsz);
        CHECK( cudaMemcpy( dmesh.elemsendudg, mesh.elemsendudg, bsz*sizeof(Int), cudaMemcpyHostToDevice ) );
    
        bsz = common.npe*common.ncAV*common.nelemsend;
        cudaTemplateMalloc(&dmesh.elemsendodg, bsz);
        CHECK( cudaMemcpy( dmesh.elemsendodg, mesh.elemsendodg, bsz*sizeof(Int), cudaMemcpyHostToDevice ) );
    }
    
    if (common.nelemrecv>0) {
        Int bsz = common.npe*common.ncu*common.nelemrecv;
        cudaTemplateMalloc(&dmesh.elemrecvind, bsz);
        CHECK( cudaMemcpy( dmesh.elemrecvind , mesh.elemrecvind , bsz*sizeof(Int), cudaMemcpyHostToDevice ) );        
    
        bsz = common.npe*common.nc*common.nelemrecv;
        cudaTemplateMalloc(&dmesh.elemrecvudg, bsz);
        CHECK( cudaMemcpy( dmesh.elemrecvudg , mesh.elemrecvudg , bsz*sizeof(Int), cudaMemcpyHostToDevice ) );
     
        bsz = common.npe*common.ncAV*common.nelemrecv;
        cudaTemplateMalloc(&dmesh.elemrecvodg, bsz);
        CHECK( cudaMemcpy( dmesh.elemrecvodg , mesh.elemrecvodg , bsz*sizeof(Int), cudaMemcpyHostToDevice ) );
    }
}

void gpuInit(solstruct &sol, resstruct &res, appstruct &app, masterstruct &master, 
       meshstruct &mesh, tempstruct &tmp, commonstruct &common, solstruct &hsol, resstruct &hres, 
       appstruct &happ, masterstruct &hmaster, meshstruct &hmesh, tempstruct &htmp, commonstruct &hcommon) 
{    
    devappstruct(app, happ);
    devmasterstruct(master, hmaster);    
    devmeshstruct(mesh, hmesh, hcommon);
    devsolstruct(sol, hsol);    
    setresstruct(res, happ, hmaster, hmesh, 0);    
    settempstruct(tmp, happ, hmaster, hmesh, 0);            
                
    // set common struct
    setcommonstruct(common, happ, hmaster, hmesh, 
            hcommon.filein, hcommon.fileout, hcommon.curvedMesh);        
    common.cpuMemory = 0;            
    
    // create cuda event handle
    CHECK(cudaEventCreate(&common.eventHandle));
    
    // create cublas handle
    CHECK_CUBLAS(cublasCreate(&common.cublasHandle));
    CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_HOST));                     //     CHECK_CUBLAS(cublasSetPointerMode(common.cublasHandle, CUBLAS_POINTER_MODE_DEVICE));    
        
    if (common.ncs>0) {
        // initialize source term
        Int N = common.npe*common.ncs*common.ne;
        cudaTemplateMalloc(&sol.sdg, N);       
        CHECK( cudaMemcpy( sol.sdg, hsol.sdg, N*sizeof(dstype), cudaMemcpyHostToDevice ) );   
        cudaTemplateMalloc(&sol.sdgg, common.nge*common.ncs*common.ne);              
    }          
    
    if (common.ncw>0) {        
        Int N = common.npe*common.ncw*common.ne;
        cudaTemplateMalloc(&sol.wsrc, N);       
        CHECK( cudaMemcpy( sol.wsrc, hsol.wsrc, N*sizeof(dstype), cudaMemcpyHostToDevice ) );    
        if (common.dae_steps>0)
            cudaTemplateMalloc(&sol.wdual, N);                 
    }          
    
    if (common.compudgavg>0) {
        cudaTemplateMalloc(&sol.udgavg, common.npe*common.nc*common.ne1+1);
    }
    
    cudaTemplateMalloc(&sol.uh, common.npf*common.ncu*common.nf);    
    #ifdef HAVE_ENZYME
        cudaTemplateMalloc(&sol.duh, common.npf*common.ncu*common.nf);
    #endif
    //cudaTemplateMalloc(&sol.uhg, common.ngf*common.ncu*common.nf);   
    //cudaTemplateMalloc(&sol.udgg, common.nge*common.nc*common.ne);   
    if (common.nco>0) {
        cudaTemplateMalloc(&sol.odgg, common.nge*common.nco*common.ne);
        cudaTemplateMalloc(&sol.og1, common.ngf*common.nco*common.nf);    
        cudaTemplateMalloc(&sol.og2, common.ngf*common.nco*common.nf);    
        #ifdef HAVE_ENZYME
        cudaTemplateMalloc(&sol.dodgg, common.nge*common.nco*common.ne);
        CHECK( cudaMemcpy( sol.dodgg, hsol.dodgg, common.nge*common.nco*common.ne*sizeof(dstype), cudaMemcpyHostToDevice ) );  

        cudaTemplateMalloc(&sol.dog1, common.ngf*common.nco*common.nf);    
        CHECK( cudaMemcpy( sol.dog1, hsol.dog1, common.ngf*common.nco*common.nf*sizeof(dstype), cudaMemcpyHostToDevice ) );  

        cudaTemplateMalloc(&sol.dog2, common.ngf*common.nco*common.nf); 
        CHECK( cudaMemcpy( sol.dog2, hsol.dog2, common.ngf*common.nco*common.nf*sizeof(dstype), cudaMemcpyHostToDevice ) );  
        #endif
    }
    
    if (common.mpiRank==0) 
        printf("finish gpuInit... \n");    
}
#endif

#endif
