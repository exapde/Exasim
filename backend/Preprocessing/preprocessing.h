#ifndef __PREPROCESSING_H__
#define __PREPROCESSING_H__

#include "structs.hpp"

class CPreprocessing {
private:
public:
    int mpirank = 0;

    ParsedSpec spec;
    InputParams params;
    PDE pde;    
    Mesh mesh;
    Master master;
    ElementClassification elemclass;
    DMD dmd;
      
    // constructor 
    CPreprocessing(string filein, int rank=0, int commsize=1); 
    
    // destructor        
    ~CPreprocessing(); 

    void SerialPreprocessing(); 

#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
    void ParallelPreprocessing(MPI_Comm comm); 
#endif  
};

#endif        

