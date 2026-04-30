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

// constructor
inline CPreprocessing::CPreprocessing(string filein, int rank, int commsize)
{
  mpirank = rank;

  params = parseInputFile(filein, rank);                           
  pde = initializePDE(params, rank);      
  pde.mpiprocs = commsize;

  spec = TextParser::parseFile(make_path(pde.datapath, pde.modelfile));        
  spec.exasimpath = pde.exasimpath;        
}

inline void CPreprocessing::SerialPreprocessing()
{  
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
    Mesh mesh = initializeParMesh(params, spec, pde, comm);   

    Master master = initializeMaster(pde, mesh, mpirank);    
    
    if (mpirank==0) {
      writepde(pde, make_path(pde.datainpath, "app.bin"));
      writemaster(master, make_path(pde.datainpath, "master.bin"));    
    }
    MPI_Barrier(comm);

    callParMetis(mesh, pde, comm); 
    
    DMD dmd = initializeDMD(mesh, master, pde, comm);     

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

