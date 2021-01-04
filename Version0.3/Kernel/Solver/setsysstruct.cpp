#ifndef __SETSYSSTRUCT
#define __SETSYSSTRUCT

void setsysstruct(sysstruct &sys, commonstruct &common, Int backend)
{
    Int ncu = common.ncu;// number of compoments of (u)    
    Int npe = common.npe; // number of nodes on master element    
    Int ne = common.ne1; // number of elements in this subdomain 
    Int N = npe*ncu*ne;    
        
    Int M;
    if (common.linearSolver==0) // GMRES
        M = common.gmresRestart+1;
    else
        M = 3;        
    M = max(M, common.RBdim);    
    
    TemplateMalloc(&sys.u, common.ndof, backend); 
    TemplateMalloc(&sys.x, common.ndof, backend); 
    TemplateMalloc(&sys.b, common.ndof, backend); 
    TemplateMalloc(&sys.r, common.ndof, backend); 
    TemplateMalloc(&sys.v, N*M, backend); 
    
    //sys.x = (dstype *) malloc(N*sizeof(dstype));
    //cout<<N<<"  "<<M<<endl;
    
    if (common.PTCparam>0) {
        TemplateMalloc(&sys.PTCmatrix, N, backend);       
        ArraySetValue(sys.PTCmatrix, one, N, backend);  
    }
    
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
}

#endif

