#ifndef __PREPROCESSING_H__
#define __PREPROCESSING_H__

#include "structs.hpp"
#include "buildstructs.hpp"

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
    CPreprocessing(std::string filein, int rank=0, int commsize=1);

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

    // HOT.7.3 — Run preprocessing into in-memory `Preprocessed`
    // bundle (no app.bin/master.bin/mesh.bin/sol.bin written). The
    // returned struct is consumed by `CSolution<M>(Preprocessed&&, ...)`.
    // Serial path (mpiprocs == 1).
    exasim::Preprocessed take();

    #if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
    // HOT.7.8 — Parallel counterpart of take(). Each rank must have
    // pre-populated `this->mesh` via meshFromArraysDistributed before
    // calling. The pipeline runs initializeMaster, buildMesh,
    // callParMetis (ParMETIS repartition), initializeDMD, then
    // builds the per-rank app/master/mesh/sol structs.
    exasim::Preprocessed takeParallel(MPI_Comm comm);
    #endif
};

#endif        

