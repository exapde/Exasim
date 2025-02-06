#ifndef __PRECONDITIONER_H__
#define __PRECONDITIONER_H__

class CPreconditioner {
private:
public:
    precondstruct precond; // store precondioner struct
        
    int mpiRank;
    
    // constructor 
    CPreconditioner(CDiscretization& disc, Int backend); 
    
    // destructor        
    ~CPreconditioner(); 

    // Construct the precontioner
    void ConstructPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend);
            
    void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int backend);
    
    void ComputeInitialGuessAndPreconditioner(sysstruct& sys, CDiscretization& disc, Int N, Int spatialScheme, Int backend);
    
    // apply the precontioner: Pv = P(u)*v
    void ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization& disc, Int backend);
    void ApplyPreconditioner(dstype* v, sysstruct& sys, CDiscretization& disc, Int spatialScheme, Int backend);
};

#endif        

