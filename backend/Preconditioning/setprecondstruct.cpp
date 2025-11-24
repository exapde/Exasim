/**
 * @file setprecondstruct.cpp
 * @brief Defines the setprecondstruct function for initializing preconditioner structures.
 *
 * @details
 * The setprecondstruct function allocates memory for the fields of a precondstruct object
 * using the specified backend. It initializes the following members:
 * - W: Workspace array of size N*M.
 * - U: Workspace array of size N*M.
 * - ipiv: Pivot array of size M+1.
 * - szipiv: Stores the size of ipiv.
 * - szW: Stores the size of W.
 * - szU: Stores the size of U.
 * - backend: Stores the backend identifier.
 *
 * @param precond Reference to the precondstruct to be initialized.
 * @param disc Reference to the CDiscretization object containing discretization information.
 * @param backend Integer specifying the memory allocation backend.
 */
#ifndef __SETPRECSTRUCT
#define __SETPRECSTRUCT

template <typename Model>
void setprecondstruct(precondstruct &precond, CDiscretization<Model>& disc, Int backend)
{    
    Int N = max(disc.common.ndof1, disc.common.ndofuhat);
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

