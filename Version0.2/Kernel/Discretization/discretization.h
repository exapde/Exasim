#ifndef __DISCRETIZATION_H__
#define __DISCRETIZATION_H__

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

    CDiscretization(){}; 
    
    // constructor for CPU 
    CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank); 
    
    // constructor for GPU 
    CDiscretization(CDiscretization& hdata); 
    
    // constructor for both CPU and GPU
    CDiscretization(string filein, string fileout, Int mpiprocs, Int mpirank, Int ompthreads, Int omprank, Int backend); 
    
    // destructor        
    ~CDiscretization(); 
        
    // compute the geometry
    void compGeometry(Int backend);    
    
    // compute the mass inverse
    void compMassInverse(Int backend);    
        
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
    
    // converge DG to CG
    void DG2CG(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
    void DG2CG2(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
    void DG2CG3(dstype* ucg, dstype* udg, dstype *utm, Int ncucg, Int ncudg, Int ncu, Int backend);
};


#endif        

