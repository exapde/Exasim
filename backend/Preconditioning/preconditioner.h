/**
 * @class CPreconditioner
 * @brief Manages the construction and application of preconditioners for linear systems.
 *
 * This class encapsulates the logic for creating and applying preconditioners, which are used to
 * improve the convergence of iterative solvers for linear systems. It interacts with discretization
 * and system structures, and supports different computational backends and spatial schemes.
 *
 * Members:
 * - precond: Stores the preconditioner structure.
 * - mpiRank: The MPI rank of the current process.
 *
 * Methods:
 * - CPreconditioner(CDiscretization& disc, Int backend): Constructor that initializes the preconditioner with the given discretization and backend.
 * - ~CPreconditioner(): Destructor.
 * - ConstructPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend): Constructs the preconditioner based on the system and discretization.
 * - ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend): Computes the initial guess and constructs the preconditioner.
 * - ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int N, Int spatialScheme, Int backend): Overloaded method to compute initial guess and preconditioner with additional parameters.
 * - ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization& disc, Int backend): Applies the preconditioner to a vector.
 * - ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization& disc, Int spatialScheme, Int backend): Overloaded method to apply the preconditioner with a specific spatial scheme.
 */
#ifndef __PRECONDITIONER_H__
#define __PRECONDITIONER_H__

template <typename Model>
class CPreconditioner {
private:
public:
    precondstruct precond; // store precondioner struct
        
    int mpiRank;
    
    // constructor 
    CPreconditioner(CDiscretization<Model>& disc, Int backend); 
    
    // destructor        
    ~CPreconditioner(); 

    // Construct the precontioner
    void ConstructPreconditioner(sysstruct& sys, CDiscretization<Model>& disc, Int backend);
            
    void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization<Model>& disc, Int backend);
    
    void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization<Model>& disc, Int N, Int spatialScheme, Int backend);
    
    // apply the precontioner: Pv = P(u)*v
    void ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization<Model>& disc, Int backend);
    void ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization<Model>& disc, Int spatialScheme, Int backend);
};

#endif        

