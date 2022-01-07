#ifndef __SOLVER_H__
#define __SOLVER_H__

class CSolver {
private:
public:
    sysstruct sys; // system struct
    
    CSolver(){}; 
    
    // constructor 
    CSolver(CDiscretization& disc, Int backend); 
    
    // destructor        
    ~CSolver(); 

    // use minimal residual to compute the initial guess
    void InitialGuess(CDiscretization& disc, CPreconditioner& prec, Int backend);
            
    // Implement PTC to solve linear/nonlinear systems
    void PseudoTransientContinuation(CDiscretization& disc, CPreconditioner& prec, Int backend);           
};

#endif        

