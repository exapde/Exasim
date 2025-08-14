#ifndef __SETSYSSTRUCT
#define __SETSYSSTRUCT

dstype rand_normal(dstype mean, dstype stddev)
{   //Box muller method
    static dstype n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        dstype x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;
            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            dstype d = sqrt(-2.0*log(r)/r);
            dstype n1 = x*d;
            n2 = y*d;
            dstype result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

void randomfield(dstype *randvect, commonstruct &common, resstruct res, meshstruct mesh, tempstruct tmp, Int backend)
{
    int N = common.npe*common.ncu*common.ne;          
    
    dstype *rvec = (dstype *) malloc((N)*sizeof(dstype));
    for (int i=0; i<N; i++) rvec[i] = rand_normal(0.0, 1.0);   
      
    //TemplateMalloc(&randvect, N, backend);   
    TemplateCopytoDevice(randvect, rvec, N, common.backend );   
            
#ifdef HAVE_MPI         
    int bsz = common.npe*common.ncu;
    
    for (int n=0; n<common.nelemsend; n++)  {       
      ArrayCopy(&tmp.tempn[bsz*n], &randvect[bsz*common.elemsend[n]], bsz);     
    }
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
    /* non-blocking send */
    Int neighbor, nsend, psend = 0, request_counter = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.tempn[psend], nsend, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    Int nrecv, precv = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.tempg[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
    
    MPI_Waitall(request_counter, common.requests, common.statuses);    
    for (int n=0; n<common.nelemrecv; n++) {        
      ArrayCopy(&randvect[bsz*common.elemrecv[n]], &tmp.tempg[bsz*n], bsz);       
    }    
#endif    
    
    Int ncu = common.ncu;
    for (Int i=0; i<ncu; i++) {
        // extract the ith component of udg and store it in res.Rq
        ArrayExtract(res.Rq, randvect, common.npe, ncu, common.ne1, 0, common.npe, i, i+1, 0, common.ne1);
        
        // make it a CG field and store in res.Ru
        ArrayDG2CG(res.Ru, res.Rq, mesh.cgent2dgent, mesh.rowent2elem, common.ndofucg);
        
        // convert CG field to DG field
        GetArrayAtIndex(res.Rq, res.Ru, mesh.cgelcon, common.npe*common.ne1);
        
        // insert utm into ucg
        ArrayInsert(randvect, res.Rq, common.npe, ncu, common.ne1, 0, common.npe, i, i+1, 0, common.ne1);
    }        
    
#ifdef HAVE_MPI             
    for (int n=0; n<common.nelemsend; n++)  {       
      ArrayCopy(&tmp.tempn[bsz*n], &randvect[bsz*common.elemsend[n]], bsz);     
    }
    
#ifdef HAVE_CUDA
    cudaDeviceSynchronize();
#endif
    
#ifdef HAVE_HIP
    hipDeviceSynchronize();
#endif
    
    /* non-blocking send */
    psend = 0;
    request_counter = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nsend = common.elemsendpts[n]*bsz;
        if (nsend>0) {
            MPI_Isend(&tmp.tempn[psend], nsend, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            psend += nsend;
            request_counter += 1;
        }
    }

    /* non-blocking receive */
    precv = 0;
    for (int n=0; n<common.nnbsd; n++) {
        neighbor = common.nbsd[n];
        nrecv = common.elemrecvpts[n]*bsz;
        if (nrecv>0) {
            MPI_Irecv(&tmp.tempg[precv], nrecv, MPI_DOUBLE, neighbor, 0,
                  MPI_COMM_WORLD, &common.requests[request_counter]);
            precv += nrecv;
            request_counter += 1;
        }
    }
    
    MPI_Waitall(request_counter, common.requests, common.statuses);    
    for (int n=0; n<common.nelemrecv; n++) {        
      ArrayCopy(&randvect[bsz*common.elemrecv[n]], &tmp.tempg[bsz*n], bsz);       
    }    
#endif        
}

void setsysstruct(sysstruct &sys, commonstruct &common, resstruct res, meshstruct mesh, tempstruct tmp, Int backend)
{
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int N = npe*ncu*ne;    
        
    Int M = common.gmresRestart+1;    
    M = max(M, common.RBdim);    
    
    // fix bug here
    Int ndof = (common.spatialScheme==0) ? N : common.ndofuhat;              
    TemplateMalloc(&sys.u, ndof, backend); 
    TemplateMalloc(&sys.x, ndof, backend); 
    TemplateMalloc(&sys.b, ndof, backend); 
    TemplateMalloc(&sys.r, ndof, backend); 
    //TemplateMalloc(&sys.v, ndof*M, backend);      
    
    if (common.spatialScheme==0) {
      TemplateMalloc(&sys.v, ndof*M, backend);      
      sys.szv = ndof * M;
    }
    else {
      sys.v = &res.K[res.szP];
      sys.szv = 0;
    }
    
    sys.backend = backend;  
    sys.szu = ndof;
    sys.szx = ndof;
    sys.szb = ndof;
    sys.szr = ndof;
    //sys.szv = ndof * M;

    ArraySetValue(sys.u, 0.0, ndof);
    ArraySetValue(sys.x, 0.0, ndof);
    ArraySetValue(sys.b, 0.0, ndof);
    ArraySetValue(sys.r, 0.0, ndof);
    ArraySetValue(sys.v, 0.0, ndof*M);
        
    if (common.ncs>0) {        
        TemplateMalloc(&sys.utmp, npe*common.nc*common.ne2, backend); 
        sys.szutmp = npe*common.nc*common.ne2;
        
        if (common.ncw>0) {
            //TemplateMalloc(&sys.w, N, backend); 
            TemplateMalloc(&sys.wtmp, npe*common.ncw*common.ne2, backend); 
            //TemplateMalloc(&sys.wsrc, N, backend);               
            sys.szwtmp = npe*common.ncw*common.ne2; 
        }                
        
        // allocate memory for the previous solutions
        if (common.temporalScheme==1) // BDF schemes 
        {
            N = common.npe*common.ncs*common.ne2;
            if (common.torder==1) {
                TemplateMalloc(&sys.udgprev1, N, backend);        
                sys.szudgprev1 = N;
            }
            else if (common.torder==2) {
                TemplateMalloc(&sys.udgprev, N, backend);      
                TemplateMalloc(&sys.udgprev1, N, backend);      
                TemplateMalloc(&sys.udgprev2, N, backend);     
                sys.szudgprev = N;
                sys.szudgprev1 = N;
                sys.szudgprev2 = N;                 
            }
            else if (common.torder==3) {
                TemplateMalloc(&sys.udgprev, N, backend);      
                TemplateMalloc(&sys.udgprev1, N, backend);      
                TemplateMalloc(&sys.udgprev2, N, backend);    
                TemplateMalloc(&sys.udgprev3, N, backend);   
                sys.szudgprev = N;
                sys.szudgprev1 = N;
                sys.szudgprev2 = N;                 
                sys.szudgprev3 = N;                  
            }      
            if (common.wave==1) {
                N = common.npe*common.ncu*common.ne1;
                if (common.torder==1) {
                    TemplateMalloc(&sys.wprev1, N, backend);   
                    sys.szwprev1 = N;
                }
                else if (common.torder==2) {
                    TemplateMalloc(&sys.wprev, N, backend);      
                    TemplateMalloc(&sys.wprev1, N, backend);      
                    TemplateMalloc(&sys.wprev2, N, backend);      
                    sys.szwprev = N;
                    sys.szwprev1 = N;
                    sys.szwprev2 = N;                
                }
                else if (common.torder==3) {
                    TemplateMalloc(&sys.wprev, N, backend);      
                    TemplateMalloc(&sys.wprev1, N, backend);      
                    TemplateMalloc(&sys.wprev2, N, backend);    
                    TemplateMalloc(&sys.wprev3, N, backend);    
                    sys.szwprev = N;
                    sys.szwprev1 = N;
                    sys.szwprev2 = N;                
                    sys.szwprev3 = N;                
                }                  
            }
        }    
        else // DIRK schemes
        {
            TemplateMalloc(&sys.udgprev, npe*common.ncs*common.ne2, backend);      
            sys.szudgprev = npe*common.ncs*common.ne2;
            if (common.ncw>0) {
                TemplateMalloc(&sys.wprev, npe*common.ncw*common.ne2, backend);                
                sys.szwprev = npe*common.ncw*common.ne2;
            }
        }        
    }    
    
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        cudaTemplateHostAlloc(&sys.tempmem, (5*M + M*M), cudaHostAllocMapped); // zero copy
        //TemplateMalloc(&sys.tempmem, (5*M + M*M), backend);            
#endif                  
    }
    else if (backend==3) { // GPU
#ifdef HAVE_HIP        
        hipTemplateHostMalloc(&sys.tempmem, (5*M + M*M), hipHostMallocMapped); // zero copy
        //TemplateMalloc(&sys.tempmem, (5*M + M*M), backend);            
#endif                  
    }    
    else { // CPU
        sys.tempmem = (dstype *) malloc((5*M + M*M)*sizeof(dstype));
    }
    sys.ipiv = (Int *) malloc(max(common.ppdegree, M*M)*sizeof(Int));             
    
    sys.szipiv = max(common.ppdegree, M*M);
    sys.sztempmem = (5*M + M*M);
        
    TemplateMalloc(&sys.randvect, common.npe*common.ncu*common.ne, backend);        
    if (common.spatialScheme==0) {
      randomfield(sys.randvect, common, res, mesh, tmp, backend);
    }
    else {
      dstype *randvectu;
      TemplateMalloc(&randvectu, common.npe*common.ncu*common.ne, backend);    
      randomfield(randvectu, common, res, mesh, tmp, backend);      
      GetFaceNodes(sys.randvect, randvectu, mesh.f2e, mesh.perm, common.npf, ncu, npe, ncu, common.nf);
    }    
    
    dstype normr = PNORM(common.cublasHandle, ndof, common.ndofuhatinterface, sys.randvect, backend);    
    //cout<<"sys.randvect: "<<common.mpiRank<<" "<<normr<<" "<<ndof<<endl;
    ArrayMultiplyScalar(common.cublasHandle, sys.randvect, 1.0/normr, ndof, backend);              
    sys.szrandvect = ndof;

    if (common.ppdegree > 1) {
        sys.lam = (dstype *) malloc((6*common.ppdegree + 2*common.ppdegree*common.ppdegree)*sizeof(dstype));        
        TemplateMalloc(&sys.q, ndof, backend);     
        TemplateMalloc(&sys.p, ndof, backend);                       
        sys.szq = ndof;
        sys.szp = ndof;
        sys.szlam = (6*common.ppdegree + 2*common.ppdegree*common.ppdegree);
    }
}

#endif

