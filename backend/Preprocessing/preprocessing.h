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
      
    // File-driven constructor: parse pdeapp.txt and pdemodel.txt.
    CPreprocessing(string filein, int rank=0, int commsize=1);

    // Programmatic constructor: skip the file parsing step entirely.
    // The caller has already populated `pde_in` / `params_in` / `spec_in`
    // (e.g. set boundary conditions, physics params, ncu/ncw, polynomial
    // order, and the Mesh inputs in InputParams). Used by external apps
    // that want to drive Exasim internals without writing a pdeapp.txt
    // to disk. The mesh itself is still fed via `pde_in.meshfile`
    // (text or binary mesh on disk) so the Mesh struct can be built
    // by initializeMesh().
    CPreprocessing(PDE pde_in, InputParams params_in, ParsedSpec spec_in,
                   int rank=0, int commsize=1);

    // destructor
    ~CPreprocessing();

    void SerialPreprocessing(); 

#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
    void ParallelPreprocessing(MPI_Comm comm); 
#endif  
};

#endif        

