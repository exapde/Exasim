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

void setsysstruct(sysstruct &sys, commonstruct &common, Int backend)
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
    TemplateMalloc(&sys.v, ndof*M, backend);         
    
    ArraySetValue(sys.u, 0.0, ndof);
    ArraySetValue(sys.x, 0.0, ndof);
    ArraySetValue(sys.b, 0.0, ndof);
    ArraySetValue(sys.r, 0.0, ndof);
    ArraySetValue(sys.v, 0.0, ndof*M);
        
    if (common.ncs>0) {        
        TemplateMalloc(&sys.utmp, npe*common.nc*common.ne2, backend); 
        
        if (common.ncw>0) {
            //TemplateMalloc(&sys.w, N, backend); 
            TemplateMalloc(&sys.wtmp, npe*common.ncw*ne, backend); 
            //TemplateMalloc(&sys.wsrc, N, backend);                         
        }                
        
        // allocate memory for the previous solutions
        if (common.temporalScheme==1) // BDF schemes 
        {
            N = common.npe*common.ncs*common.ne2;
            if (common.torder==1) {
                TemplateMalloc(&sys.udgprev1, N, backend);                    
            }
            else if (common.torder==2) {
                TemplateMalloc(&sys.udgprev, N, backend);      
                TemplateMalloc(&sys.udgprev1, N, backend);      
                TemplateMalloc(&sys.udgprev2, N, backend);                      
            }
            else if (common.torder==3) {
                TemplateMalloc(&sys.udgprev, N, backend);      
                TemplateMalloc(&sys.udgprev1, N, backend);      
                TemplateMalloc(&sys.udgprev2, N, backend);    
                TemplateMalloc(&sys.udgprev3, N, backend);    
            }      
            if (common.wave==1) {
                N = common.npe*common.ncu*common.ne1;
                if (common.torder==1) {
                    TemplateMalloc(&sys.wprev1, N, backend);   
                }
                else if (common.torder==2) {
                    TemplateMalloc(&sys.wprev, N, backend);      
                    TemplateMalloc(&sys.wprev1, N, backend);      
                    TemplateMalloc(&sys.wprev2, N, backend);                      
                }
                else if (common.torder==3) {
                    TemplateMalloc(&sys.wprev, N, backend);      
                    TemplateMalloc(&sys.wprev1, N, backend);      
                    TemplateMalloc(&sys.wprev2, N, backend);    
                    TemplateMalloc(&sys.wprev3, N, backend);    
                }                  
            }
        }    
        else // DIRK schemes
        {
            TemplateMalloc(&sys.udgprev, npe*common.ncs*common.ne2, backend);      
            if (common.ncw>0) 
                TemplateMalloc(&sys.wprev, npe*common.ncw*ne, backend);                
        }        
    }    
    
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        cudaTemplateHostAlloc(&sys.tempmem, (5*M + M*M), cudaHostAllocMapped); // zero copy
        //TemplateMalloc(&sys.tempmem, (5*M + M*M), backend);    
#endif        
        sys.cpuMemory = 0;    
    }
    else { // CPU
        sys.tempmem = (dstype *) malloc((5*M + M*M)*sizeof(dstype));
        sys.cpuMemory = 1;    
    }
    sys.ipiv = (Int *) malloc(max(common.ppdegree, M*M)*sizeof(Int));             
    
    // sys.normcu = (dstype *) malloc(ncu*sizeof(dstype));    
        
    N = ndof; // fix bug here                
    TemplateMalloc(&sys.randvect, N, backend);    
#ifdef HAVE_CUDA                               
    dstype *rvec = (dstype *) malloc((N)*sizeof(dstype));
    for (int i=0; i<N; i++) rvec[i] = rand_normal(0.0, 1.0);   
    CHECK( cudaMemcpy(sys.randvect, rvec, N*sizeof(dstype), cudaMemcpyHostToDevice ) );  
    free(rvec);
#endif                
#ifndef HAVE_CUDA      
    for (int i=0; i<N; i++) sys.randvect[i] = rand_normal(0.0, 1.0);        
#endif                   
    dstype normr = PNORM(common.cublasHandle, N, sys.randvect, backend);    
//    cout<<common.mpiRank<<" "<<normr<<" "<<N<<endl;
    ArrayMultiplyScalar(common.cublasHandle, sys.randvect, 1.0/normr, N, backend);                       
//     normr = PNORM(common.cublasHandle, N, sys.randvect, backend);    
//     cout<<common.mpiRank<<" "<<normr<<" "<<N<<endl;
//     error("here");
    
    if (common.ppdegree > 1) {
        sys.lam = (dstype *) malloc((6*common.ppdegree + 2*common.ppdegree*common.ppdegree)*sizeof(dstype));        
        TemplateMalloc(&sys.q, ndof, backend);     
        TemplateMalloc(&sys.p, ndof, backend);                       
    }
}

#endif

