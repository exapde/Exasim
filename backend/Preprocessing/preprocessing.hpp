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

// HOT.7.3 — programmatic preprocessing entry point.
//
// Runs the same pipeline as `SerialPreprocessing` (initializeMesh +
// initializeMaster + buildMesh + initializeDMD + apply_bcm) but
// returns a typed `Preprocessed` bundle instead of writing
// app.bin/master.bin/mesh.bin/sol.bin to disk. The caller hands
// the bundle to `CSolution<M>(Preprocessed&&, ...)`.
//
// Currently serial-only (mpiprocs == 1). MPI variant lands in 7.6.
inline exasim::Preprocessed CPreprocessing::take()
{
    using exasim::Preprocessed;

    // Mesh + master setup (file-driven path runs these inside
    // SerialPreprocessing; the programmatic path can pre-populate
    // `mesh` via meshFromArrays before the call).
    if (mesh.np == 0)
        mesh = initializeMesh(params, pde);
    Master mas = initializeMaster(pde, mesh, mpirank);

    // Mirrors writeBinaryFiles' spec → pde adjustment so dimensions
    // are consistent with what the legacy file path computes.
    for (const auto& vec : spec.vectors) {
        const std::string& name = vec.first;
        int size = vec.second;
        if (name == "uhat") pde.ncu = size;
        if (name == "v")    pde.ncv = size;
        if (name == "w")    pde.ncw = size;
        if (name == "uq")   pde.nc  = size;
    }
    for (size_t i = 0; i < spec.functions.size(); ++i) {
        const auto& fn = spec.functions[i];
        if (fn.name == "VisScalars") pde.nsca  = fn.outputsize;
        if (fn.name == "VisVectors") pde.nvec  = fn.outputsize / pde.nd;
        if (fn.name == "VisTensors") pde.nten  = fn.outputsize / (pde.nd*pde.nd);
        if (fn.name == "QoIboundary") pde.nsurf = fn.outputsize;
        if (fn.name == "QoIvolume")  pde.nvqoi = fn.outputsize;
    }

    if (pde.mpiprocs > 1)
        error("CPreprocessing::take(): MPI path lands in HOT.7.6. Use SerialPreprocessing or pde.mpiprocs=1.");

    // buildMesh populates mesh.xdg, mesh.f, mesh.t2lf, mesh.localfaces.
    buildMesh(mesh, pde, mas);

    // Serial DMD + bf — same as the writeBinaryFiles serial branch.
    DMD dmd_local = initializeDMD(pde, mesh);
    int  ne_local = (int)dmd_local.elempart.size();
    std::vector<int> bf((size_t)mesh.nfe * ne_local);
    {
        std::vector<int> fi((size_t)mesh.nfe * ne_local);
        select_columns(fi.data(), mesh.f.data(), dmd_local.elempart.data(), mesh.nfe, ne_local);
        apply_bcm(bf.data(), fi.data(), mesh.boundaryConditions.data(),
                  mesh.nfe*ne_local, mesh.nbcm);
    }

    // Build the four runtime structs directly from typed sources.
    Preprocessed out;
    out.app    = exasim::buildAppStruct(pde);
    out.master = exasim::buildMasterStruct(mas);
    out.mesh   = exasim::buildMeshStruct(mesh, mas, dmd_local, bf,
                                         pde.hybrid, pde.mpiprocs, out.ti);
    out.sol    = exasim::buildSolStruct(mesh, mas, pde, dmd_local);

    // Free the char-array boundary expression buffers — they were
    // moved into mesh by initializeMesh and are no longer needed
    // (the Mesh struct's destructor doesn't own them).
    if (mesh.nbndexpr > 0) freeCharArray(mesh.boundaryExprs, mesh.nbndexpr);
    if (mesh.nbndexpr > 0) freeCharArray(mesh.curvedBoundaryExprs, mesh.nbndexpr);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs1, mesh.nprdexpr*mesh.nprdcom);
    if (mesh.nprdexpr > 0) freeCharArray(mesh.periodicExprs2, mesh.nprdexpr*mesh.nprdcom);

    return out;
}

#endif        

