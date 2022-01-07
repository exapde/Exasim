#ifndef __SOLVER
#define __SOLVER

#include "solver.h"
#include "setsysstruct.cpp"
#include "conjgrad.cpp"
#include "gmres.cpp"
#include "lmsolver.cpp"
#include "ptcsolver.cpp"

// constructor
CSolver::CSolver(CDiscretization& disc, Int backend)
{
    // set up sys struct
    setsysstruct(sys, disc.common, backend);
    
    // copy wdg to w for wave problem
    //if (disc.common.wave==1)
    //    ArrayCopy(sys.w, disc.sol.wdg, disc.common.ndof1, backend);    
}

// destructor
CSolver::~CSolver()
{    
    sys.freememory(sys.cpuMemory);    
}

void CSolver::InitialGuess(CDiscretization& disc, CPreconditioner& prec, Int backend)
{   
    LevenbergMarquardt(sys, disc, prec, backend); 
}

void CSolver::PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, Int backend)       
{
    // solve the system using PTC to obtain the solution sys.u
    disc.common.nonlinearSolverIter = PTCsolver(sys, disc, prec, backend); 
    
    // update UDG
    disc.updateUDG(sys.u, backend);           
}

#endif        

