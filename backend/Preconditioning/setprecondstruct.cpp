#ifndef __SETPRECSTRUCT
#define __SETPRECSTRUCT

void setprecondstruct(precondstruct &precond, CDiscretization& disc, Int backend)
{    
    Int N = disc.common.ndof1;
    Int M = disc.common.RBdim;    
    
    TemplateMalloc(&precond.W, N*M, backend); 
    TemplateMalloc(&precond.U, N*M, backend);   
    //TemplateMalloc(&precond.V, N*M, backend);   
    //TemplateMalloc(&precond.R, N*M, backend);   
    //TemplateMalloc(&precond.H, M*M, backend);          
    //TemplateMalloc(&precond.y, N, backend);     
    //precond.y = &disc.res.Ru[0];   
    TemplateMalloc(&precond.ipiv, M+1, backend);                
    
    precond.szipiv = M+1;
    precond.szW = N*M;
    precond.szU = N*M;

    // if (disc.common.preconditioner>0)
    //     TemplateMalloc(&precond.z, N, backend); 
            
    // if (disc.common.precMatrixType==0) {
    //     TemplateMalloc(&precond.C, N, backend); 
    //     ArraySetValue(precond.C, one, N);  
    // }
    
    precond.backend = backend;       
}

#endif

