#ifndef __SETPRECSTRUCT
#define __SETPRECSTRUCT

void setprecondstruct(precondstruct &precond, CDiscretization& disc, Int backend)
{    
    Int N = disc.common.ndof1;
    Int M = disc.common.RBdim;    
    
    TemplateMalloc(&precond.W, N*M, backend); 
    TemplateMalloc(&precond.U, N*M, backend);   
    TemplateMalloc(&precond.ipiv, M+1, backend);                
    
    precond.szipiv = M+1;
    precond.szW = N*M;
    precond.szU = N*M;
    
    precond.backend = backend;       
}

#endif

