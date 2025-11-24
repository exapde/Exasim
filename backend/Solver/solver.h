/**
 * @class CSolver
 * @brief Solver class for handling linear and nonlinear systems in the Exasim backend.
 *
 * This class encapsulates the solver logic, including pseudo-transient continuation (PTC)
 * and Newton-based methods for solving systems arising from discretization.
 *
 * Members:
 * - sys: System structure containing problem-specific data.
 * - mpiRank: MPI rank for parallel computations.
 *
 * Constructors & Destructors:
 * - CSolver(CDiscretization& disc, Int backend): Constructs a solver with given discretization and backend.
 * - ~CSolver(): Destructor for cleanup.
 *
 * Methods:
 * - void PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int backend):
 *     Implements the PTC method to solve linear or nonlinear systems.
 * - void NewtonSolver(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int N, Int spatialScheme, Int backend):
 *     Solves systems using the Newton method with configurable parameters.
 */
#ifndef __SOLVER_H__
#define __SOLVER_H__

template <typename Model>
class CSolver {
private:
public:
    sysstruct sys; // system struct
    
    int mpiRank;
    
    // constructor 
    CSolver(CDiscretization<Model>& disc, Int backend); 
    
    // destructor        
    ~CSolver(); 
            
    // Implement PTC to solve linear/nonlinear systems
    void PseudoTransientContinuation(CDiscretization<Model>& disc, CPreconditioner<Model>& prec, ofstream& out, Int backend);           

    void NewtonSolver(CDiscretization<Model>& disc, CPreconditioner<Model>& prec, ofstream& out, Int N, Int spatialScheme, Int backend);       
};

#endif        

