#ifndef __PREPROCESSING
#define __PREPROCESSING

#ifdef HAVE_PARMETIS
  #include <metis.h>
#ifdef HAVE_MPI
  #include <parmetis.h>
#endif
#else
  #include <cstdint>
  using idx_t = std::int32_t;
#endif

#include "preprocessing.h"

#include "TextParser.hpp"
#include "tinyexpr.hpp"
#include "helpersexasim.hpp"
#include "readpdeapp.hpp"
#include "readmesh.hpp"
#include "makemeshexasim.hpp"
#include "makemaster.hpp"
#include "makemasterexasim.hpp"
#include "domaindecomposition.hpp"
#include "writebinaryfilesexasim.hpp"

#ifdef HAVE_PARMETIS
#ifdef HAVE_MPI
#include "parmetisexasim.hpp"
#endif
#endif

// File-driven constructor.
inline CPreprocessing::CPreprocessing(string filein, int rank, int commsize)
{
  mpirank = rank;

  params = parseInputFile(filein, rank);
  pde = initializePDE(params, rank);
  pde.mpiprocs = commsize;

  spec = TextParser::parseFile(make_path(pde.datapath, pde.modelfile));
  spec.exasimpath = pde.exasimpath;
}

// Programmatic constructor: caller passes pre-populated structs and we
// skip parseInputFile / initializePDE / TextParser::parseFile entirely.
// We still call pdeFinalizeDerived() so the runtime-side flag /
// problem / factor / solversparam vectors are packed (otherwise
// readsolstruct dereferences a NULL flag pointer).
inline CPreprocessing::CPreprocessing(PDE pde_in, InputParams params_in,
                                      ParsedSpec spec_in,
                                      int rank, int commsize)
{
  mpirank = rank;
  params  = std::move(params_in);
  pde     = std::move(pde_in);
  pde.mpiprocs = commsize;
  spec    = std::move(spec_in);
  if (spec.exasimpath.empty()) spec.exasimpath = pde.exasimpath;
  pdeFinalizeDerived(pde);
}

inline void CPreprocessing::SerialPreprocessing()
{
    // If `mesh` was injected from outside (e.g. via meshFromArrays /
    // setMesh) skip the file-read step. Otherwise read mesh from disk
    // per `pde.meshfile`. The legacy file-driven path leaves `mesh.np
    // == 0` after the constructor and falls through to initializeMesh.
    if (mesh.np == 0)
        mesh = initializeMesh(params, pde);

    master = initializeMaster(pde, mesh);
    writeBinaryFiles(pde, mesh, master, spec);

    if (mesh.nbndexpr > 0) freeCharArray(mesh.boundaryExprs, mesh.nbndexpr);
    if (mesh.nbndexpr > 0) freeCharArray(mesh.curvedBoundaryExprs, mesh.nbndexpr);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs1, mesh.nprdexpr*mesh.nprdcom);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs2, mesh.nprdexpr*mesh.nprdcom);
}

#if defined(HAVE_PARMETIS) && defined(HAVE_MPI)
inline void CPreprocessing::ParallelPreprocessing(MPI_Comm comm)
{
    // Same mesh-injection convention as SerialPreprocessing: if the
    // caller pre-populated `this->mesh` (programmatic path), skip the
    // file read. The parallel-mesh-from-arrays path lands in a later
    // slice — for now the legacy file path runs.
    if (mesh.np == 0)
        mesh = initializeParMesh(params, spec, pde, comm);

    master = initializeMaster(pde, mesh, mpirank);

    if (mpirank==0) {
      writepde(pde, make_path(pde.datainpath, "app.bin"));
      writemaster(master, make_path(pde.datainpath, "master.bin"));
    }
    MPI_Barrier(comm);

    callParMetis(mesh, pde, comm);

    dmd = initializeDMD(mesh, master, pde, comm);

    writemesh(mesh, dmd, pde, master, comm);

    writesol(mesh, dmd, pde, master, comm);

    if (mesh.nbndexpr > 0) freeCharArray(mesh.boundaryExprs, mesh.nbndexpr);
    if (mesh.nbndexpr > 0) freeCharArray(mesh.curvedBoundaryExprs, mesh.nbndexpr);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs1, mesh.nprdexpr*mesh.nprdcom);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs2, mesh.nprdexpr*mesh.nprdcom);
}
#endif

// destructor
inline CPreprocessing::~CPreprocessing()
{            
    if (mpirank==0) printf("CPreprocessing destructor is called successfully.\n");
}

#endif        

