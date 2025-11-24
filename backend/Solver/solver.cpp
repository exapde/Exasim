/*
    CSolver class

    This class provides methods for solving nonlinear systems arising from discretized PDEs.
    It interfaces with system structures, discretization, and preconditioners, and supports
    different solution strategies.

    Members:
    - CSolver(CDiscretization& disc, Int backend)
        Constructor. Initializes the solver system structure using the provided discretization
        and backend. Optionally prints system information if in debug mode.

    - ~CSolver()
        Destructor. Frees memory allocated for the system structure and prints a message
        if running on the root MPI rank.

    - void PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int backend)
        Solves the nonlinear system using the Pseudo-Transient Continuation (PTC) method.
        Updates the solution vector and the discretization's UDG field.

    - void NewtonSolver(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int N, Int spatialScheme, Int backend)
        Solves the nonlinear system using a Newton-based solver. The number of iterations (N)
        and the spatial discretization scheme can be specified. Updates the UDG field if
        spatialScheme is zero.

    Dependencies:
    - solver.h
    - setsysstruct.cpp
    - getpoly.cpp
    - gmres.cpp
    - ptcsolver.cpp

    Notes:
    - Uses MPI for parallelism; certain actions are performed only on the root rank.
    - Assumes existence of system structure 'sys' and appropriate member functions.
*/
#ifndef __SOLVER
#define __SOLVER

#include "solver.h"
#include "setsysstruct.cpp"
#include "getpoly.cpp"
#include "gmres.cpp"
#include "ptcsolver.cpp"

// constructor
template <typename Model>
CSolver<Model>::CSolver(CDiscretization<Model>& disc, Int backend)
{
    mpiRank = disc.common.mpiRank;
    setsysstruct(sys, disc.common, disc.res, disc.mesh, disc.tmp, backend);   
    if ((disc.common.mpiRank==0) && (disc.common.debugMode==1)) sys.printinfo(); 

    if (disc.common.mpiRank == 0) printf("finish CSolver constructor... \n");    
}

// destructor
template <typename Model>
CSolver<Model>::~CSolver()
{        
    sys.freememory(sys.backend);    
    if (mpiRank==0) printf("CSolver destructor: sys memory is freed successfully.\n");
}

template <typename Model>
void CSolver<Model>::PseudoTransientContinuation(CDiscretization<Model>& disc, CPreconditioner<Model>& prec, ofstream& out, Int backend)       
{
    // solve the system using PTC to obtain the solution sys.u
    disc.common.nonlinearSolverIter = PTCsolver(sys, disc, prec, out, backend); 
    
    // update UDG
    disc.updateUDG(sys.u, backend);           
}

template <typename Model>
void CSolver<Model>::NewtonSolver(CDiscretization<Model>& disc, CPreconditioner<Model>& prec, ofstream& out, Int N, Int spatialScheme, Int backend)       
{
    // solve the system using PTC to obtain the solution sys.u
    disc.common.nonlinearSolverIter = NonlinearSolver(sys, disc, prec, out, N, spatialScheme, backend); 
    
    // update UDG
    if (spatialScheme==0) disc.updateUDG(sys.u, backend);           
}

#endif        

