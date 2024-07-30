#ifndef __SOLVER_H__
#define __SOLVER_H__

class CSolver {
private:
public:
    sysstruct sys; // system struct
    
    int mpiRank;
    
    // constructor 
    CSolver(CDiscretization& disc, Int backend); 
    
    // destructor        
    ~CSolver(); 
            
    // Implement PTC to solve linear/nonlinear systems
    void PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int backend);           

    void NewtonSolver(CDiscretization& disc, CPreconditioner& prec, ofstream& out, Int N, Int spatialScheme, Int backend);       
};

#endif        

