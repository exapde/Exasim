/**
 * @class CDiscretization
 * @brief Handles the discretization process for numerical simulations, supporting both CPU and GPU backends.
 *
 * This class encapsulates data structures and methods required for discretizing PDEs, assembling linear systems,
 * evaluating residuals, fluxes, artificial viscosity fields, and performing conversions between DG and CG representations.
 *
 * Members:
 * - sol: Solution structure containing state variables.
 * - res: Residual structure for storing residuals.
 * - app: Application-specific parameters and data.
 * - master: Master element data for discretization.
 * - mesh: Mesh structure containing grid information.
 * - tmp: Temporary storage for intermediate computations.
 * - common: Common parameters shared across computations.
 *
 * Methods:
 * - CDiscretization(...): Constructor initializing the discretization with input/output files and parallelization parameters.
 * - ~CDiscretization(): Destructor for cleanup.
 * - compGeometry(...): Computes geometry-related quantities.
 * - compMassInverse(...): Computes the inverse of the mass matrix.
 * - hdgAssembleLinearSystem(...): Assembles the linear system for HDG methods.
 * - hdgAssembleResidual(...): Assembles the residual vector for HDG methods.
 * - evalResidual(...): Evaluates the residual vector.
 * - evalResidual(...): Evaluates the residual vector at a given solution.
 * - evalQ(...): Evaluates the flux and stores it in the solution structure.
 * - evalQSer(...): Serial version of flux evaluation.
 * - evalQ(...): Evaluates the flux at a given solution.
 * - evalMatVec(...): Evaluates the matrix-vector product Jv = J(u)*v.
 * - evalMatVec(...): Overloaded version with spatial scheme selection.
 * - updateUDG(...): Updates the solution structure with new values and computes flux.
 * - updateU(...): Updates the solution structure with new values.
 * - evalAVfield(...): Evaluates the artificial viscosity field at a given solution.
 * - evalAVfield(...): Evaluates the artificial viscosity field at the current solution.
 * - evalOutput(...): Evaluates output quantities at the current solution.
 * - evalMonitor(...): Evaluates a monitor function for tracking solution changes during pseudotime stepping.
 * - DG2CG(...): Converts DG representation to CG.
 * - DG2CG2(...): Alternative DG to CG conversion.
 * - DG2CG3(...): Another variant of DG to CG conversion.
 */
#ifndef __DISCRETIZATION_H__
#define __DISCRETIZATION_H__


template <typename Model>
class CDiscretization {
private:
public:
    solstruct sol;
    resstruct res;
    appstruct app;
    masterstruct master; 
    meshstruct mesh;
    tempstruct tmp;    
    commonstruct common;
    // solstruct hsol;

    // constructor for both CPU and GPU
    CDiscretization(string filein, string fileout, string exasimpath, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank, Int backend); 
    
    // destructor        
    ~CDiscretization(); 
        
    // compute the geometry
    void compGeometry(Int backend);    
    
    // compute the mass inverse
    void compMassInverse(Int backend);    

    void hdgAssembleLinearSystem(dstype *b, Int backend);        
    void hdgAssembleResidual(dstype *b, Int backend);        

    // evaluate the residual vector 
    void evalResidual(Int backend);
        
    // evaluate the residual vector at u
    void evalResidual(dstype* Ru, dstype* u, Int backend);
    
    // evaluate the flux q and store it in sol.udg
    void evalQ(Int backend);
    
    // evaluate the flux q and store it in sol.udg
    void evalQSer(Int backend);
    
    // evaluate the flux q at u
    void evalQ(dstype* q, dstype* u, Int backend);
        
    // evaluate the matrix-vector product Jv = J(u)*v
    void evalMatVec(dstype* Jv, dstype* v, dstype* u, dstype* Ru, Int backend);   

    void evalMatVec(dstype* Jv, dstype* v, dstype* u, dstype* Ru, Int spatialScheme, Int backend);   
    
    // insert u into sol.udg and compute q
    void updateUDG(dstype* u, Int backend);

    // insert u into sol.udg
    void updateU(dstype* u, Int backend);    
    
    // evaluate the artificial viscosity field at u
    void evalAVfield(dstype* avField, dstype* u, Int backend);       
        
    // evaluate the artificial viscosity field at sol.udg
    void evalAVfield(dstype* avField, Int backend);       
    
    // evaluate the artificial viscosity field at sol.udg
    void evalOutput(dstype* output, Int backend);       
    
    // evaluate a monitor function to monitor changes in solution QoIs for pseudotime stepping
    void evalMonitor(dstype* output, dstype* udg, dstype* wdg, Int nc, Int backed);
    
    // converge DG to CG
    void DG2CG(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
    void DG2CG2(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
    void DG2CG3(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
};


#endif        

