#ifndef __SOLVER
#define __SOLVER

#include "solver.h"
#include "setsysstruct.cpp"
#include "getpoly.cpp"
#include "gmres.cpp"
#include "ptcsolver.cpp"

// constructor
CSolver::CSolver(CDiscretization& disc, Int backend)
{
    mpiRank = disc.common.mpiRank;
    setsysstruct(sys, disc.common, disc.res, disc.mesh, disc.tmp, backend);   
    if ((disc.common.mpiRank==0) && (disc.common.debugMode==1)) sys.printinfo(); 
}

// destructor
CSolver::~CSolver()
{        
    sys.freememory(sys.backend);    
    if (mpiRank==0) printf("CSolver destructor: sys memory is freed successfully.\n");
}

void CSolver::PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int backend)       
{
    // solve the system using PTC to obtain the solution sys.u
    disc.common.nonlinearSolverIter = PTCsolver(sys, disc, prec, out, backend); 
    
    // update UDG
    disc.updateUDG(sys.u, backend);           
}

void CSolver::NewtonSolver(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int N, Int spatialScheme, Int backend)       
{
    // solve the system using PTC to obtain the solution sys.u
    disc.common.nonlinearSolverIter = NonlinearSolver(sys, disc, prec, out, N, spatialScheme, backend); 
    
    // update UDG
    if (spatialScheme==0) disc.updateUDG(sys.u, backend);           
}

#endif        

