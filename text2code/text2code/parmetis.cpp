#include <mpi.h>
#include <parmetis.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

/// Parallel mesh partitioning using ParMETIS_V3_PartMeshKway
///
/// elmdist: global element distribution over MPI ranks, size = nprocs + 1.
/// e2n:    local element-to-node connectivity (flattened, global node IDs)
///
/// Output:
///   epart_local[e] = partition ID of local element e
///
void partitionMeshParMETIS(std::vector<idx_t>       &epart_local,
                           std::vector<idx_t>       &e2n,      // local connectivity
                           const std::vector<idx_t> &elmdist,
                           idx_t                     nve,
                           idx_t                     ncommon,
                           idx_t                     nparts,
                           MPI_Comm                  comm)
{
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    //----------------------------------------------------------------------
    // 1. Local number of elements
    //----------------------------------------------------------------------
    idx_t ne_local = elmdist[rank+1] - elmdist[rank];

    if (ne_local > 0 && e2n.size() != (size_t)ne_local * nve) {
        if (rank == 0) {
            std::cerr << "partitionMeshParMETIS: e2n.size() = " << e2n.size()
                      << " but expected ne_local*nve = " << (ne_local * nve)
                      << std::endl;
        }
        return;
    }

    //----------------------------------------------------------------------
    // 2. Build eptr = element pointer array
    //----------------------------------------------------------------------
    std::vector<idx_t> eptr(ne_local + 1);
    eptr[0] = 0;
    for (idx_t e = 0; e < ne_local; ++e)
        eptr[e+1] = eptr[e] + nve;

    //----------------------------------------------------------------------
    // 3. Allocate epart_local directly as ParMETIS output buffer
    //----------------------------------------------------------------------
    epart_local.resize(ne_local);    // ParMETIS writes partition IDs here

    //----------------------------------------------------------------------
    // 4. ParMETIS parameters
    //----------------------------------------------------------------------
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t ncon    = 1;
    idx_t edgecut = 0;
    idx_t options[3] = {0, 0, 0};

    std::vector<real_t> tpwgts(nparts * ncon, real_t(1.0) / real_t(nparts));
    std::vector<real_t> ubvec(ncon, real_t(1.05));

    //----------------------------------------------------------------------
    // 5. Raw pointers required by the C API
    //----------------------------------------------------------------------
    idx_t *elmdist_ptr = const_cast<idx_t*>(elmdist.data());
    idx_t *eptr_ptr    = eptr.data();
    idx_t *eind_ptr    = e2n.data();

    idx_t *part_ptr    = epart_local.data();   // *** Direct output buffer ***

    //----------------------------------------------------------------------
    // 6. Call ParMETIS
    //----------------------------------------------------------------------
    int status = ParMETIS_V3_PartMeshKway(
        elmdist_ptr,
        eptr_ptr,
        eind_ptr,
        nullptr,           // elmwgt
        &wgtflag,
        &numflag,
        &ncon,
        &ncommon,
        &nparts,
        tpwgts.data(),
        ubvec.data(),
        options,
        &edgecut,
        part_ptr,          // wrote directly into epart_local
        &comm
    );

    if (status != METIS_OK) {
        if (rank == 0)
            std::cerr << "ParMETIS_V3_PartMeshKway failed with status "
                      << status << std::endl;
        return;
    }

    if (rank == 0) {
        std::cout << "Finished partitioning mesh using ParMETIS (edgecut = "
                  << edgecut << ")" << std::endl;
    }
}

std::vector<idx_t> buildElmdistBalanced(idx_t ne_global, MPI_Comm comm)
{
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    std::vector<idx_t> elmdist(size + 1);
    elmdist[0] = 0;

    idx_t base = ne_global / size;
    idx_t rem  = ne_global % size;

    for (int p = 0; p < size; ++p) {
        idx_t count = base + (p < rem ? 1 : 0);
        elmdist[p+1] = elmdist[p] + count;
    }

    // Everyone should see the same elmdist (we built it deterministically).
    return elmdist;
}

// Assumes idx_t is the ParMETIS/METIS index type (usually 32-bit or 64-bit integer).
// If idx_t is 64-bit, replace MPI_INT below with MPI_LONG_LONG accordingly.
std::vector<idx_t> buildElmdistFromLocalCount(idx_t ne_local, MPI_Comm comm)
{
    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    // Gather ne_local from all ranks
    std::vector<idx_t> ne_per_rank(size);

    // If idx_t is 64-bit, use MPI_LONG_LONG here instead of MPI_INT
    MPI_Allgather(&ne_local, 1, MPI_INT,
                  ne_per_rank.data(), 1, MPI_INT,
                  comm);

    // Build elmdist: prefix sum of ne_per_rank
    std::vector<idx_t> elmdist(size + 1);
    elmdist[0] = 0;
    for (int p = 0; p < size; ++p) {
        elmdist[p+1] = elmdist[p] + ne_per_rank[p];
    }

    // Optional sanity check: sum over ranks
    idx_t ne_global = elmdist[size];
    if (rank == 0) {
        std::cout << "buildElmdistFromLocalCount: ne_global = "
                  << ne_global << std::endl;
    }

    return elmdist;
}

// Simple contiguous partition helper
inline void computeLocalRange(int globalN, int size, int rank,
                              int &localN, int &offset)
{
    int base = globalN / size;
    int rem  = globalN % size;
    if (rank < rem) {
        localN = base + 1;
        offset = rank * localN;
    } else {
        localN = base;
        offset = rem * (base + 1) + (rank - rem) * base;
    }
}


// -----------------------------------------------------------------------------
// Parallel reader: each rank reads *its portion* of p and t from the same file
// -----------------------------------------------------------------------------
void readParMeshFromBinaryFile(const std::string& filename,
                               Mesh& mesh,
                               MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // -----------------------------------------------------------------
    // 1. Rank 0 reads the global header (nd, np_global, nve, ne_global)
    // -----------------------------------------------------------------
    int nd_global = 0, np_global = 0, nve_global = 0, ne_global = 0;

    if (rank == 0) {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);

        int* ndims = readiarrayfromdouble(in, 4);
        nd_global  = ndims[0];
        np_global  = ndims[1];
        nve_global = ndims[2];
        ne_global  = ndims[3];
        free(ndims);

        in.close();
    }

    // Broadcast global mesh sizes to all ranks
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // -----------------------------------------------------------------
    // 2. Compute local np and ne and their offsets in the global arrays
    // -----------------------------------------------------------------
    int np_local   = 0;
    int nodeOffset = 0;  // starting node index (0-based) for this rank

    int ne_local   = 0;
    int elemOffset = 0;  // starting element index (0-based) for this rank

    computeLocalRange(np_global, size, rank, np_local, nodeOffset);
    computeLocalRange(ne_global, size, rank, ne_local, elemOffset);

    // -----------------------------------------------------------------
    // 3. Compute byte offsets into the file
    //
    // File layout (as implied by your serial reader):
    //   - 4 doubles: nd, np, nve, ne
    //   - p: np_global * nd_global doubles
    //   - t: ne_global * nve_global doubles (each storing an integer)
    // -----------------------------------------------------------------
    const std::size_t headerDoubles = 4;
    const std::size_t pDoubles      = static_cast<std::size_t>(np_global) * nd_global;
    const std::size_t tDoubles      = static_cast<std::size_t>(ne_global) * nve_global;

    const std::size_t headerBytes = headerDoubles * sizeof(double);
    const std::size_t pBytes      = pDoubles      * sizeof(double);
    const std::size_t tBytes      = tDoubles      * sizeof(double);

    // Our local slices:
    const std::size_t localPDoubles = static_cast<std::size_t>(np_local) * nd_global;
    const std::size_t localTDoubles = static_cast<std::size_t>(ne_local) * nve_global;

    const std::size_t localPBytes = localPDoubles * sizeof(double);
    const std::size_t localTBytes = localTDoubles * sizeof(double);

    // Start offsets (in doubles) for this rank
    const std::size_t nodeStartDouble = static_cast<std::size_t>(nodeOffset) * nd_global;
    const std::size_t elemStartDouble = static_cast<std::size_t>(elemOffset) * nve_global;

    const std::size_t pStartByte = headerBytes + nodeStartDouble * sizeof(double);
    const std::size_t tStartByte = headerBytes + pBytes + elemStartDouble * sizeof(double);

    // -----------------------------------------------------------------
    // 4. Resize local mesh storage and fill metadata
    // -----------------------------------------------------------------
    mesh.nd  = nd_global;   // same on all ranks
    mesh.nve = nve_global;  // same on all ranks
    mesh.np  = np_local;    // local
    mesh.ne  = ne_local;    // local

    mesh.p.resize(static_cast<std::size_t>(mesh.np) * mesh.nd);
    mesh.t.resize(static_cast<std::size_t>(mesh.ne) * mesh.nve);

    if (np_local == 0 && ne_local == 0) {
        // This rank owns no data; nothing more to do
        return;
    }

    // -----------------------------------------------------------------
    // 5. Each rank opens the file and reads its slices of p and t
    // -----------------------------------------------------------------
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) error("Unable to open file " + filename);

    // 5a. Read local p (coordinates)
    if (np_local > 0) {
        in.seekg(static_cast<std::streamoff>(pStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(mesh.p.data()),
                static_cast<std::streamsize>(localPBytes));
        if (!in) error("Error reading local p block from " + filename);
    }

    // 5b. Read local t (connectivity), stored as doubles in file
    if (ne_local > 0) {
        in.seekg(static_cast<std::streamoff>(tStartByte), std::ios::beg);

        std::vector<double> tbuf(localTDoubles);
        in.read(reinterpret_cast<char*>(tbuf.data()),
                static_cast<std::streamsize>(localTBytes));
        if (!in) error("Error reading local t block from " + filename);

        // Convert double -> int, consistent with readiarrayfromdouble
        for (std::size_t i = 0; i < localTDoubles; ++i) {
            mesh.t[i] = static_cast<int>(tbuf[i]);
        }
    }

    in.close();
}

// ----------------------------------------------------------
// Generic parallel mesh reader:
//  - ReaderFunc must be: void(const std::string&, Mesh&)
//  - tagBase: MPI tag base (tagBase for p, tagBase+1 for t)
// ----------------------------------------------------------
template <typename ReaderFunc>
void readParMeshGeneric(const std::string& filename,
                        Mesh& mesh,
                        MPI_Comm comm,
                        ReaderFunc reader,
                        int tagBase)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    Mesh globalMesh;   // only used on rank 0

    int nd_global  = 0;
    int np_global  = 0;
    int nve_global = 0;
    int ne_global  = 0;

    // 1. Rank 0 reads full mesh using provided reader
    if (rank == 0) {
        reader(filename, globalMesh);
        nd_global  = globalMesh.nd;
        np_global  = globalMesh.np;
        nve_global = globalMesh.nve;
        ne_global  = globalMesh.ne;
    }

    // Broadcast global metadata
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // 2. Compute local np, ne and offsets for this rank
    int np_local   = 0;
    int nodeOffset = 0;
    int ne_local   = 0;
    int elemOffset = 0;

    computeLocalRange(np_global, size, rank, np_local, nodeOffset);
    computeLocalRange(ne_global, size, rank, ne_local, elemOffset);

    // 3. Initialize local mesh metadata and allocate storage
    mesh.nd  = nd_global;
    mesh.nve = nve_global;
    mesh.np  = np_local;
    mesh.ne  = ne_local;

    mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
    mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);

    if (np_local == 0 && ne_local == 0) {
        // This rank owns nothing; done
        return;
    }

    // 4. Distribute nodal coordinates (mesh.p), column-major: p(j + i*nd)
    const int tagP = tagBase;
    if (rank == 0) {
        // Rank 0: fill its local block directly
        for (int i = 0; i < np_local; ++i) {
            int gnode = nodeOffset + i;
            for (int j = 0; j < nd_global; ++j) {
                mesh.p[j + i * nd_global] =
                    globalMesh.p[j + gnode * nd_global];
            }
        }

        // Send other ranks their blocks
        for (int r = 1; r < size; ++r) {
            int np_r, nodeOff_r;
            computeLocalRange(np_global, size, r, np_r, nodeOff_r);
            if (np_r == 0) continue;

            std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
            for (int i = 0; i < np_r; ++i) {
                int gnode = nodeOff_r + i;
                for (int j = 0; j < nd_global; ++j) {
                    sendP[j + i * nd_global] =
                        globalMesh.p[j + gnode * nd_global];
                }
            }

            MPI_Send(sendP.data(),
                     static_cast<int>(sendP.size()),
                     MPI_DOUBLE, r, tagP, comm);
        }
    } else {
        // Non-root ranks receive their p block
        if (np_local > 0) {
            std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
            MPI_Recv(recvP.data(),
                     static_cast<int>(recvP.size()),
                     MPI_DOUBLE, 0, tagP, comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < np_local; ++i)
                for (int j = 0; j < nd_global; ++j)
                    mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
        }
    }

    // 5. Distribute connectivity (mesh.t), column-major: t(j + i*nve)
    const int tagT = tagBase + 1;
    if (rank == 0) {
        // Rank 0: fill its local block
        for (int i = 0; i < ne_local; ++i) {
            int gele = elemOffset + i;
            for (int j = 0; j < nve_global; ++j) {
                mesh.t[j + i * nve_global] =
                    globalMesh.t[j + gele * nve_global];
            }
        }

        // Send other ranks their blocks
        for (int r = 1; r < size; ++r) {
            int ne_r, elemOff_r;
            computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
            if (ne_r == 0) continue;

            std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
            for (int i = 0; i < ne_r; ++i) {
                int gele = elemOff_r + i;
                for (int j = 0; j < nve_global; ++j) {
                    sendT[j + i * nve_global] =
                        globalMesh.t[j + gele * nve_global];
                }
            }

            MPI_Send(sendT.data(),
                     static_cast<int>(sendT.size()),
                     MPI_INT, r, tagT, comm);
        }
    } else {
        // Non-root ranks receive their t block
        if (ne_local > 0) {
            std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
            MPI_Recv(recvT.data(),
                     static_cast<int>(recvT.size()),
                     MPI_INT, 0, tagT, comm, MPI_STATUS_IGNORE);

            for (int i = 0; i < ne_local; ++i)
                for (int j = 0; j < nve_global; ++j)
                    mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
        }
    }
}

// Text
void readParMeshFromTextFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 500, 501
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromTextFile,
                       500);
}

// VTK
void readParMeshFromVTKFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 300, 301
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromVTKFile,
                       300);
}

// VTU
void readParMeshFromVTUFile(const std::string& filename,
                            Mesh& mesh,
                            MPI_Comm comm)
{
    // tags 400, 401
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromVTUFile,
                       400);
}

// Gmsh v2
void readParMeshFromMshV2File(const std::string& filename,
                              Mesh& mesh,
                              MPI_Comm comm)
{
    // tags 100, 101
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromMshV2File,
                       100);
}

// Gmsh v4
void readParMeshFromMshV4File(const std::string& filename,
                              Mesh& mesh,
                              MPI_Comm comm)
{
    // tags 200, 201
    readParMeshGeneric(filename, mesh, comm,
                       readMeshFromMshV4File,
                       200);
}

void readParMeshFromGmshFile(const std::string& filename, Mesh& mesh, MPI_Comm comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    double version = 0.0;
    int fileType   = 0;
    int foundMeshFormat = 0; // bool as int for MPI

    if (rank == 0) {
        std::ifstream infile(filename);
        if (!infile.is_open()) error("Error opening file: " + filename);

        std::string line;
        while (std::getline(infile, line)) {
            if (line == "$MeshFormat") {
                foundMeshFormat = 1;
                break;
            }
        }

        if (!foundMeshFormat) {
            std::cerr << "Error in readParMeshFromGmshFile: $MeshFormat not found in file: "
                      << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::getline(infile, line); // read version line
        std::istringstream versionStream(line);
        versionStream >> version >> fileType;

        infile.close();
    }

    // Broadcast detection results to all ranks
    MPI_Bcast(&foundMeshFormat, 1, MPI_INT, 0, comm);
    MPI_Bcast(&version,         1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&fileType,        1, MPI_INT,    0, comm);

    if (!foundMeshFormat) {
        // Rank 0 already printed an error; make sure everyone stops.
        std::cerr << "Error in readParMeshFromGmshFile: $MeshFormat not found (broadcast)\n";
        std::exit(EXIT_FAILURE);
    }

    if (fileType != 0) {
        if (rank == 0) {
            std::cerr << "Error in readParMeshFromGmshFile: Only ASCII .msh files are supported (binary detected)"
                      << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }

    if (version >= 4.0) {
        if (rank == 0)
            std::cout << "Detected Gmsh format v" << version << " → using parallel V4 reader\n";
        readParMeshFromMshV4File(filename, mesh, comm);
    } else if (version >= 2.0 && version < 4.0) {
        if (rank == 0)
            std::cout << "Detected Gmsh format v" << version << " → using parallel V2/V3 reader\n";
        readParMeshFromMshV2File(filename, mesh, comm);
    } else {
        if (rank == 0) {
            std::cerr << "Error in readParMeshFromGmshFile: Unsupported Gmsh mesh format version "
                      << version << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }
}

void readParMeshFromFile(const std::string& filename, Mesh& mesh, MPI_Comm comm) 
{
    std::string ext = getFileExtension(filename);

    if (ext == "bin") {
        readMeshFromBinaryFile(filename, mesh, comm);
    } else if (ext == "txt") {
        readMeshFromTextFile(filename, mesh, comm);
    } else if (ext == "msh") {
        readMeshFromGmshFile(filename, mesh, comm);
    } else if (ext == "vtk") {
        readMeshFromVTKFile(filename, mesh, comm);
    } else if (ext == "vtu") {
        readMeshFromVTUFile(filename, mesh, comm);
    } else {
        error("Unsupported mesh file format: " + ext);        
    }
        
    for (int i=0; i<mesh.ne * mesh.nve; i++) mesh.t[i] -= 1;      
    
    // determine element type
    if (mesh.nd == 2 && mesh.nve == 4) mesh.elemtype = 1;
    else if (mesh.nd == 3 && mesh.nve == 8) mesh.elemtype = 1;
    else mesh.elemtype = 0;

    mesh.nvf = (mesh.nd == 3) ? (mesh.nd + mesh.elemtype) : mesh.nd;
    mesh.nfe = mesh.nd + (mesh.nd - 1) * mesh.elemtype + 1;                    
}

void readLocalMeshFromBinaryFile(const std::string& filename,
                                 Mesh& mesh,
                                 const std::vector<idx_t>& epart_local,
                                 MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // -----------------------------------------------------------------
    // 1. Rank 0 reads global header (nd, np_global, nve, ne_global)
    // -----------------------------------------------------------------
    int nd_global = 0, np_global = 0, nve_global = 0, ne_global = 0;

    if (rank == 0) {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in) error("Unable to open file " + filename);

        int* ndims = readiarrayfromdouble(in, 4);
        nd_global  = ndims[0];
        np_global  = ndims[1];
        nve_global = ndims[2];
        ne_global  = ndims[3];
        free(ndims);

        in.close();
    }

    // Broadcast global mesh sizes to all ranks
    int header[4] = { nd_global, np_global, nve_global, ne_global };
    MPI_Bcast(header, 4, MPI_INT, 0, comm);
    nd_global  = header[0];
    np_global  = header[1];
    nve_global = header[2];
    ne_global  = header[3];

    // -----------------------------------------------------------------
    // 2. Reconstruct original block distribution used for ParMETIS
    //    (same computeLocalRange as in readParMeshFromBinaryFile)
    // -----------------------------------------------------------------
    int ne_local_orig = 0;
    int elemOffset    = 0;  // starting global element index for this rank
    computeLocalRange(ne_global, size, rank, ne_local_orig, elemOffset);

    // Sanity check: local ParMETIS output must match original local ne
    if (static_cast<int>(epart_local.size()) != ne_local_orig) {
        std::cerr << "readLocalMeshFromBinaryFile: epart_local.size() = "
                  << epart_local.size()
                  << " but original ne_local = " << ne_local_orig
                  << " on rank " << rank << "\n";
        error("Inconsistent local ParMETIS partition size");
    }

    // -----------------------------------------------------------------
    // 3. Build global element-partition array epart_global[e] (size ne_global)
    // -----------------------------------------------------------------
    std::vector<idx_t> epart_global(ne_global);

    // Determine MPI datatype for idx_t
    MPI_Datatype mpi_idx_t = (sizeof(idx_t) == sizeof(int))
                             ? MPI_INT
                             : MPI_LONG_LONG; // adjust if you use 64-bit idx_t differently

    // Allgatherv parameters
    std::vector<int> recvcounts(size), displs(size);
    for (int r = 0; r < size; ++r) {
        int ne_r = 0, off_r = 0;
        computeLocalRange(ne_global, size, r, ne_r, off_r);
        recvcounts[r] = ne_r;
        displs[r]     = off_r;
    }

    MPI_Allgatherv(epart_local.data(),
                   ne_local_orig,
                   mpi_idx_t,
                   epart_global.data(),
                   recvcounts.data(),
                   displs.data(),
                   mpi_idx_t,
                   comm);

    // -----------------------------------------------------------------
    // 4. Determine which global elements belong to *this* rank after ParMETIS
    // -----------------------------------------------------------------
    std::vector<idx_t> ownedElems; // global element indices for this rank
    ownedElems.reserve(ne_global / size + 1);

    for (idx_t e = 0; e < (idx_t)ne_global; ++e) {
        if (epart_global[e] == rank) {
            ownedElems.push_back(e);
        }
    }

    const int ne_local = static_cast<int>(ownedElems.size());

    // If this rank owns no elements, also owns no nodes
    if (ne_local == 0) {
        mesh.nd  = nd_global;
        mesh.nve = nve_global;
        mesh.ne  = 0;
        mesh.np  = 0;
        mesh.p.clear();
        mesh.t.clear();
        return;
    }

    // -----------------------------------------------------------------
    // 5. Compute file layout constants
    //
    // File layout:
    //   - 4 doubles: nd, np, nve, ne
    //   - p: np_global * nd_global doubles
    //   - t: ne_global * nve_global doubles (each storing an integer)
    // -----------------------------------------------------------------
    const std::size_t headerDoubles = 4;
    const std::size_t pDoubles      = static_cast<std::size_t>(np_global) * nd_global;
    const std::size_t tDoubles      = static_cast<std::size_t>(ne_global) * nve_global;

    const std::size_t headerBytes   = headerDoubles * sizeof(double);
    const std::size_t pBytes        = pDoubles      * sizeof(double);
    (void)tDoubles;

    const std::size_t tBlockStartByte = headerBytes + pBytes;

    // -----------------------------------------------------------------
    // 6. Read connectivity for owned elements and build node mask
    // -----------------------------------------------------------------
    mesh.nd  = nd_global;
    mesh.nve = nve_global;
    mesh.ne  = ne_local;

    // Temporarily store global node IDs
    mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);

    std::vector<char> nodeMask(np_global, 0);

    std::ifstream in(filename, std::ios::in | std::ios::binary);
    if (!in) error("Unable to open file " + filename);

    std::vector<double> tmpConn(nve_global);

    for (int eL = 0; eL < ne_local; ++eL) {
        idx_t eGlobal = ownedElems[eL];

        std::size_t elemStartDouble = static_cast<std::size_t>(eGlobal) * nve_global;
        std::size_t elemStartByte   = tBlockStartByte + elemStartDouble * sizeof(double);

        in.seekg(static_cast<std::streamoff>(elemStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(tmpConn.data()),
                static_cast<std::streamsize>(nve_global * sizeof(double)));
        if (!in) error("Error reading connectivity for element from " + filename);

        for (int j = 0; j < nve_global; ++j) {
            int gnode = static_cast<int>(tmpConn[j]);
            mesh.t[j + eL * nve_global] = gnode; // global node for now

            if (gnode < 0 || gnode >= np_global) {
                error("readLocalMeshFromBinaryFile: node index out of range");
            }
            nodeMask[gnode] = 1;
        }
    }

    // -----------------------------------------------------------------
    // 7. Build list of local nodes and global→local mapping
    // -----------------------------------------------------------------
    std::vector<int> localNodes;
    localNodes.reserve(np_global);

    std::vector<int> g2l(np_global, -1);

    for (int g = 0; g < np_global; ++g) {
        if (nodeMask[g]) {
            int l = static_cast<int>(localNodes.size());
            localNodes.push_back(g);
            g2l[g] = l;
        }
    }

    const int np_local = static_cast<int>(localNodes.size());
    mesh.np = np_local;

    // -----------------------------------------------------------------
    // 8. Read coordinates for local nodes
    // -----------------------------------------------------------------
    mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);

    // Clear EOF flags before re-using the stream
    in.clear();

    std::vector<double> coordBuf(nd_global);

    for (int l = 0; l < np_local; ++l) {
        int gnode = localNodes[l];

        std::size_t nodeStartDouble = static_cast<std::size_t>(gnode) * nd_global;
        std::size_t nodeStartByte   = headerBytes + nodeStartDouble * sizeof(double);

        in.seekg(static_cast<std::streamoff>(nodeStartByte), std::ios::beg);
        in.read(reinterpret_cast<char*>(coordBuf.data()),
                static_cast<std::streamsize>(nd_global * sizeof(double)));
        if (!in) error("Error reading coordinates for node from " + filename);

        for (int j = 0; j < nd_global; ++j) {
            mesh.p[j + l * nd_global] = coordBuf[j];
        }
    }

    in.close();

    // -----------------------------------------------------------------
    // 9. Remap connectivity from global node IDs to local node IDs
    // -----------------------------------------------------------------
    for (int eL = 0; eL < ne_local; ++eL) {
        for (int j = 0; j < nve_global; ++j) {
            int gnode = mesh.t[j + eL * nve_global];
            int lnode = g2l[gnode];
            if (lnode < 0) {
                error("readLocalMeshFromBinaryFile: missing local mapping for node");
            }
            mesh.t[j + eL * nve_global] = lnode;
        }
    }
}

#include <mpi.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>

// Your Mesh struct
struct Mesh {
    std::vector<double> p; // nd × np (column-major)
    std::vector<int>    t; // nve × ne (column-major)
    int nd;   // spatial dimension
    int np;   // local number of nodes
    int nve;  // vertices per element
    int ne;   // local number of elements
};

void error(const std::string& msg);

// Helper: compute prefix offsets from counts
inline void prefixSums(const std::vector<int>& counts, std::vector<int>& displs, int& total)
{
    displs.resize(counts.size());
    total = 0;
    for (std::size_t i = 0; i < counts.size(); ++i) {
        displs[i] = total;
        total += counts[i];
    }
}

// Main migration routine
void migrateMeshWithParMETIS(const Mesh& mesh_in,
                             const std::vector<idx_t>& epart_local,
                             Mesh& mesh_out,
                             MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int nd  = mesh_in.nd;
    const int nve = mesh_in.nve;
    const int ne_local_in = mesh_in.ne;
    const int np_local_in = mesh_in.np;

    if ((int)epart_local.size() != ne_local_in) {
        error("migrateMeshWithParMETIS: epart_local size does not match local ne");
    }

    // ------------------------------------------------------------------
    // 1. ELEMENT MIGRATION (based on ParMETIS partition epart_local)
    // ------------------------------------------------------------------

    // 1a. Count how many elements we send to each destination rank
    std::vector<int> sendElemCounts(size, 0);
    for (int e = 0; e < ne_local_in; ++e) {
        int dest = static_cast<int>(epart_local[e]);
        if (dest < 0 || dest >= size) {
            error("migrateMeshWithParMETIS: epart_local has invalid destination rank");
        }
        sendElemCounts[dest] += 1;
    }

    // 1b. Exchange counts to know how many elements we will receive
    std::vector<int> recvElemCounts(size, 0);
    MPI_Alltoall(sendElemCounts.data(), 1, MPI_INT,
                 recvElemCounts.data(), 1, MPI_INT,
                 comm);

    // 1c. Build displacements and total send/recv element counts
    std::vector<int> sendElemDispls, recvElemDispls;
    int totalSendElems = 0, totalRecvElems = 0;
    prefixSums(sendElemCounts, sendElemDispls, totalSendElems);
    prefixSums(recvElemCounts, recvElemDispls, totalRecvElems);

    // 1d. Pack connectivity to send: global node IDs, nve per element
    std::vector<int> sendElemConn(static_cast<std::size_t>(totalSendElems) * nve);

    // per-destination "cursor" in element space
    std::vector<int> elemOffsetPerDest(size, 0);

    for (int e = 0; e < ne_local_in; ++e) {
        int dest = static_cast<int>(epart_local[e]);
        int localIdxInDest = elemOffsetPerDest[dest]++;
        int globalElemPos  = sendElemDispls[dest] + localIdxInDest; // element index in send buffer

        // Copy connectivity column e → position globalElemPos
        for (int j = 0; j < nve; ++j) {
            int gnode = mesh_in.t[j + e * nve]; // global node ID
            sendElemConn[j + globalElemPos * nve] = gnode;
        }
    }

    // 1e. Receive connectivity for new local elements
    std::vector<int> recvElemConn(static_cast<std::size_t>(totalRecvElems) * nve);

    // For Alltoallv, counts are in *number of ints*, not elements
    std::vector<int> sendElemCountsInts(size), recvElemCountsInts(size);
    for (int r = 0; r < size; ++r) {
        sendElemCountsInts[r] = sendElemCounts[r] * nve;
        recvElemCountsInts[r] = recvElemCounts[r] * nve;
        sendElemDispls[r]    *= nve;
        recvElemDispls[r]    *= nve;
    }

    MPI_Alltoallv(sendElemConn.data(),
                  sendElemCountsInts.data(), sendElemDispls.data(), MPI_INT,
                  recvElemConn.data(),
                  recvElemCountsInts.data(), recvElemDispls.data(), MPI_INT,
                  comm);

    // Now recvElemConn holds all elements that this rank will own (global node IDs)
    const int ne_local_out = totalRecvElems;

    // ------------------------------------------------------------------
    // 2. BUILD NODE OWNERSHIP INFORMATION (from initial distribution)
    //    We assume nodes were initially distributed contiguously:
    //      rank r owns global nodes in [nodedist[r], nodedist[r+1])
    // ------------------------------------------------------------------
    // Gather np_local_in from all ranks to reconstruct nodedist
    std::vector<int> np_per_rank(size);
    MPI_Allgather(&np_local_in, 1, MPI_INT,
                  np_per_rank.data(), 1, MPI_INT,
                  comm);

    std::vector<int> nodedist(size + 1);
    nodedist[0] = 0;
    for (int r = 0; r < size; ++r) {
        nodedist[r + 1] = nodedist[r] + np_per_rank[r];
    }
    const int np_global = nodedist[size];

    // Helper to find owner rank of global node g
    auto findNodeOwner = [&](int gnode) -> int {
        // simple linear search over ranks; size is usually modest
        for (int r = 0; r < size; ++r) {
            if (gnode >= nodedist[r] && gnode < nodedist[r + 1]) {
                return r;
            }
        }
        return -1; // should never happen
    };

    // ------------------------------------------------------------------
    // 3. DETERMINE WHICH GLOBAL NODES THIS RANK NEEDS FOR ITS NEW ELEMENTS
    // ------------------------------------------------------------------
    // 3a. Mark all global nodes referenced by new local elements
    std::vector<char> nodeMask(np_global, 0);

    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            int gnode = recvElemConn[j + e * nve];
            if (gnode < 0 || gnode >= np_global) {
                error("migrateMeshWithParMETIS: global node index out of range");
            }
            nodeMask[gnode] = 1;
        }
    }

    // 3b. For each owner rank, collect the global nodes we need from it
    std::vector<std::vector<int>> needNodes(size); // needNodes[owner] = list of GIDs
    for (int g = 0; g < np_global; ++g) {
        if (nodeMask[g]) {
            int owner = findNodeOwner(g);
            if (owner < 0) {
                error("migrateMeshWithParMETIS: could not find owner for global node");
            }
            needNodes[owner].push_back(g);
        }
    }

    // Total #local nodes after migration
    int np_local_out = 0;
    for (int r = 0; r < size; ++r) {
        np_local_out += static_cast<int>(needNodes[r].size());
    }

    // ------------------------------------------------------------------
    // 4. REQUEST NODE COORDINATES FROM OWNERS (two-phase MPI exchange)
    // ------------------------------------------------------------------

    // 4a. Phase 1: send requests for node IDs
    std::vector<int> sendNodeReqCounts(size, 0);
    for (int r = 0; r < size; ++r) {
        sendNodeReqCounts[r] = static_cast<int>(needNodes[r].size());
    }

    std::vector<int> recvNodeReqCounts(size, 0);
    MPI_Alltoall(sendNodeReqCounts.data(), 1, MPI_INT,
                 recvNodeReqCounts.data(), 1, MPI_INT,
                 comm);

    std::vector<int> sendNodeReqDispls, recvNodeReqDispls;
    int totalSendNodeReq = 0, totalRecvNodeReq = 0;
    prefixSums(sendNodeReqCounts, sendNodeReqDispls, totalSendNodeReq);
    prefixSums(recvNodeReqCounts, recvNodeReqDispls, totalRecvNodeReq);

    // Flatten our node requests into a buffer
    std::vector<int> sendNodeReqBuf(totalSendNodeReq);
    for (int r = 0; r < size; ++r) {
        int offset = sendNodeReqDispls[r];
        for (std::size_t k = 0; k < needNodes[r].size(); ++k) {
            sendNodeReqBuf[offset + static_cast<int>(k)] = needNodes[r][k];
        }
    }

    // Owners receive node ID requests
    std::vector<int> recvNodeReqBuf(totalRecvNodeReq);
    MPI_Alltoallv(sendNodeReqBuf.data(),
                  sendNodeReqCounts.data(), sendNodeReqDispls.data(), MPI_INT,
                  recvNodeReqBuf.data(),
                  recvNodeReqCounts.data(), recvNodeReqDispls.data(), MPI_INT,
                  comm);

    // 4b. Phase 2: owners respond with coordinates in the same order

    // Each owner will send: recvNodeReqCounts[dest] * nd doubles
    std::vector<int> sendCoordCounts(size, 0);
    for (int r = 0; r < size; ++r) {
        sendCoordCounts[r] = recvNodeReqCounts[r] * nd;
    }

    std::vector<int> recvCoordCounts(size, 0);
    MPI_Alltoall(sendCoordCounts.data(), 1, MPI_INT,
                 recvCoordCounts.data(), 1, MPI_INT,
                 comm);

    std::vector<int> sendCoordDispls, recvCoordDispls;
    int totalSendCoord = 0, totalRecvCoord = 0;
    prefixSums(sendCoordCounts, sendCoordDispls, totalSendCoord);
    prefixSums(recvCoordCounts, recvCoordDispls, totalRecvCoord);

    std::vector<double> sendCoordBuf(totalSendCoord);

    // Fill sendCoordBuf: for each destination, in the order we received requests
    for (int src = 0; src < size; ++src) {
        int countNodesFromSrc = recvNodeReqCounts[src]; // nodes requested by 'src' from this rank
        int reqOffset         = recvNodeReqDispls[src];
        int coordOffset       = sendCoordDispls[src];

        for (int k = 0; k < countNodesFromSrc; ++k) {
            int gnode = recvNodeReqBuf[reqOffset + k];

            // Map global node ID to local index on this owner
            int ownerLocalStart = nodedist[rank];
            int ownerLocalEnd   = nodedist[rank + 1];

            if (gnode < ownerLocalStart || gnode >= ownerLocalEnd) {
                error("migrateMeshWithParMETIS: node request to wrong owner");
            }

            int lidx = gnode - ownerLocalStart; // local index in mesh_in.p

            // Copy nd coordinates
            for (int j = 0; j < nd; ++j) {
                sendCoordBuf[coordOffset + k * nd + j] =
                    mesh_in.p[j + lidx * nd];
            }
        }
    }

    // Perform coordinate exchange
    std::vector<double> recvCoordBuf(totalRecvCoord);
    MPI_Alltoallv(sendCoordBuf.data(),
                  sendCoordCounts.data(), sendCoordDispls.data(), MPI_DOUBLE,
                  recvCoordBuf.data(),
                  recvCoordCounts.data(), recvCoordDispls.data(), MPI_DOUBLE,
                  comm);

    // ------------------------------------------------------------------
    // 5. BUILD FINAL LOCAL NODE LIST & GLOBAL→LOCAL MAPPING
    // ------------------------------------------------------------------

    mesh_out.nd  = nd;
    mesh_out.nve = nve;
    mesh_out.ne  = ne_local_out;
    mesh_out.np  = np_local_out;

    mesh_out.t.resize(static_cast<std::size_t>(nve) * ne_local_out);
    mesh_out.p.resize(static_cast<std::size_t>(nd)  * np_local_out);

    // 5a. Build owner-based offsets into local node space
    std::vector<int> ownerNodeOffsets(size);
    {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            ownerNodeOffsets[r] = offset;
            offset += static_cast<int>(needNodes[r].size());
        }
    }

    // 5b. Global→local index map
    std::vector<int> g2l(np_global, -1);
    for (int r = 0; r < size; ++r) {
        int base = ownerNodeOffsets[r];
        for (std::size_t k = 0; k < needNodes[r].size(); ++k) {
            int gnode = needNodes[r][k];
            int lnode = base + static_cast<int>(k);
            g2l[gnode] = lnode;
        }
    }

    // 5c. Place received coordinates into mesh_out.p
    for (int owner = 0; owner < size; ++owner) {
        int numNodesFromOwner = sendNodeReqCounts[owner]; // nodes we requested from 'owner'
        if (numNodesFromOwner == 0) continue;

        int coordOffset = recvCoordDispls[owner]; // in doubles
        int baseLocal   = ownerNodeOffsets[owner];

        for (int k = 0; k < numNodesFromOwner; ++k) {
            int localNode = baseLocal + k;
            for (int j = 0; j < nd; ++j) {
                mesh_out.p[j + localNode * nd] =
                    recvCoordBuf[coordOffset + k * nd + j];
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. REMAP CONNECTIVITY FROM GLOBAL NODE IDS → LOCAL NODE IDS
    // ------------------------------------------------------------------

    // First copy new element connectivity (still global IDs) into mesh_out.t
    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            mesh_out.t[j + e * nve] = recvElemConn[j + e * nve];
        }
    }

    // Then replace global IDs with local IDs using g2l
    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            int gnode = mesh_out.t[j + e * nve];
            int lnode = g2l[gnode];
            if (lnode < 0) {
                error("migrateMeshWithParMETIS: missing local mapping for node");
            }
            mesh_out.t[j + e * nve] = lnode;
        }
    }
}

void migrateMeshWithParMETIS(const Mesh& mesh_in,
                             const std::vector<idx_t>& epart_local,
                             Mesh& mesh_out,
                             MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int nd  = mesh_in.nd;
    const int nve = mesh_in.nve;
    const int ne_local_in = mesh_in.ne;
    const int np_local_in = mesh_in.np;

    if ((int)epart_local.size() != ne_local_in) {
        error("migrateMeshWithParMETIS: epart_local size does not match local ne");
    }

    // ------------------------------------------------------------------
    // 0. RECONSTRUCT GLOBAL ELEMENT IDS FOR CURRENT DISTRIBUTION
    //    We assume initial elements are distributed contiguously:
    //      rank r owns global elements [elmdist[r], elmdist[r+1])
    // ------------------------------------------------------------------
    // Gather ne_local_in from all ranks
    std::vector<int> ne_per_rank(size);
    MPI_Allgather(&ne_local_in, 1, MPI_INT,
                  ne_per_rank.data(), 1, MPI_INT,
                  comm);

    std::vector<int> elmdist(size + 1);
    elmdist[0] = 0;
    for (int r = 0; r < size; ++r) {
        elmdist[r + 1] = elmdist[r] + ne_per_rank[r];
    }
    const int ne_global = elmdist[size];

    // Global ID of local element e on this rank:
    const int gstart = elmdist[rank];
    std::vector<int> elemGlobalID_in(ne_local_in);
    for (int e = 0; e < ne_local_in; ++e) {
        elemGlobalID_in[e] = gstart + e;
    }

    // ------------------------------------------------------------------
    // 1. ELEMENT MIGRATION (based on ParMETIS partition epart_local)
    // ------------------------------------------------------------------

    // 1a. Count how many elements we send to each destination rank
    std::vector<int> sendElemCounts(size, 0);
    for (int e = 0; e < ne_local_in; ++e) {
        int dest = static_cast<int>(epart_local[e]);
        if (dest < 0 || dest >= size) {
            error("migrateMeshWithParMETIS: epart_local has invalid destination rank");
        }
        sendElemCounts[dest] += 1;
    }

    // 1b. Exchange counts to know how many elements we will receive
    std::vector<int> recvElemCounts(size, 0);
    MPI_Alltoall(sendElemCounts.data(), 1, MPI_INT,
                 recvElemCounts.data(), 1, MPI_INT,
                 comm);

    // 1c. Build displacements and total send/recv element counts
    std::vector<int> sendElemDispls, recvElemDispls;
    int totalSendElems = 0, totalRecvElems = 0;
    prefixSums(sendElemCounts, sendElemDispls, totalSendElems);
    prefixSums(recvElemCounts, recvElemDispls, totalRecvElems);

    // 1d. Pack connectivity and global element IDs to send
    std::vector<int> sendElemConn(static_cast<std::size_t>(totalSendElems) * nve);
    std::vector<int> sendElemGlobalID(static_cast<std::size_t>(totalSendElems));

    // per-destination "cursor" in element space
    std::vector<int> elemOffsetPerDest(size, 0);

    for (int e = 0; e < ne_local_in; ++e) {
        int dest           = static_cast<int>(epart_local[e]);
        int localIdxInDest = elemOffsetPerDest[dest]++;
        int elemPos        = sendElemDispls[dest] + localIdxInDest; // element index in send buffers

        // Copy connectivity column e → position elemPos
        for (int j = 0; j < nve; ++j) {
            int gnode = mesh_in.t[j + e * nve]; // global node ID
            sendElemConn[j + elemPos * nve] = gnode;
        }

        // Copy global element ID
        sendElemGlobalID[elemPos] = elemGlobalID_in[e];
    }

    // 1e. Receive connectivity for new local elements
    std::vector<int> recvElemConn(static_cast<std::size_t>(totalRecvElems) * nve);
    std::vector<int> recvElemGlobalID(static_cast<std::size_t>(totalRecvElems));

    // For Alltoallv, counts are in *number of ints*, not elements
    std::vector<int> sendElemCountsInts(size), recvElemCountsInts(size);
    for (int r = 0; r < size; ++r) {
        sendElemCountsInts[r] = sendElemCounts[r] * nve;
        recvElemCountsInts[r] = recvElemCounts[r] * nve;
        sendElemDispls[r]    *= nve;
        recvElemDispls[r]    *= nve;
    }

    // Connectivity exchange
    MPI_Alltoallv(sendElemConn.data(),
                  sendElemCountsInts.data(), sendElemDispls.data(), MPI_INT,
                  recvElemConn.data(),
                  recvElemCountsInts.data(), recvElemDispls.data(), MPI_INT,
                  comm);

    // Global element ID exchange (one int per element)
    // Use original element counts/displs (in units of "elements")
    MPI_Alltoallv(sendElemGlobalID.data(),
                  sendElemCounts.data(),   // counts in elements
                  sendElemDispls.data(),   // displs in elements
                  MPI_INT,
                  recvElemGlobalID.data(),
                  recvElemCounts.data(),
                  recvElemDispls.data(),
                  MPI_INT,
                  comm);

    // Now recvElemConn holds all elements this rank will own (global node IDs),
    // and recvElemGlobalID[e] is the global element ID of new local element e.
    const int ne_local_out = totalRecvElems;

    // ------------------------------------------------------------------
    // 2. BUILD NODE OWNERSHIP INFORMATION (from initial distribution)
    //    We assume nodes were initially distributed contiguously:
    //      rank r owns global nodes in [nodedist[r], nodedist[r+1])
    // ------------------------------------------------------------------
    // Gather np_local_in from all ranks to reconstruct nodedist
    std::vector<int> np_per_rank(size);
    MPI_Allgather(&np_local_in, 1, MPI_INT,
                  np_per_rank.data(), 1, MPI_INT,
                  comm);

    std::vector<int> nodedist(size + 1);
    nodedist[0] = 0;
    for (int r = 0; r < size; ++r) {
        nodedist[r + 1] = nodedist[r] + np_per_rank[r];
    }
    const int np_global = nodedist[size];

    // Helper to find owner rank of global node g
    auto findNodeOwner = [&](int gnode) -> int {
        for (int r = 0; r < size; ++r) {
            if (gnode >= nodedist[r] && gnode < nodedist[r + 1]) {
                return r;
            }
        }
        return -1; // should never happen
    };

    // ------------------------------------------------------------------
    // 3. DETERMINE WHICH GLOBAL NODES THIS RANK NEEDS FOR ITS NEW ELEMENTS
    // ------------------------------------------------------------------
    std::vector<char> nodeMask(np_global, 0);

    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            int gnode = recvElemConn[j + e * nve];
            if (gnode < 0 || gnode >= np_global) {
                error("migrateMeshWithParMETIS: global node index out of range");
            }
            nodeMask[gnode] = 1;
        }
    }

    // For each owner rank, collect the global nodes we need from it
    std::vector<std::vector<int>> needNodes(size); // needNodes[owner] = list of GIDs
    for (int g = 0; g < np_global; ++g) {
        if (nodeMask[g]) {
            int owner = findNodeOwner(g);
            if (owner < 0) {
                error("migrateMeshWithParMETIS: could not find owner for global node");
            }
            needNodes[owner].push_back(g);
        }
    }

    // Total #local nodes after migration
    int np_local_out = 0;
    for (int r = 0; r < size; ++r) {
        np_local_out += static_cast<int>(needNodes[r].size());
    }

    // ------------------------------------------------------------------
    // 4. REQUEST NODE COORDINATES FROM OWNERS (two-phase MPI exchange)
    // ------------------------------------------------------------------

    // 4a. Phase 1: send requests for node IDs
    std::vector<int> sendNodeReqCounts(size, 0);
    for (int r = 0; r < size; ++r) {
        sendNodeReqCounts[r] = static_cast<int>(needNodes[r].size());
    }

    std::vector<int> recvNodeReqCounts(size, 0);
    MPI_Alltoall(sendNodeReqCounts.data(), 1, MPI_INT,
                 recvNodeReqCounts.data(), 1, MPI_INT,
                 comm);

    std::vector<int> sendNodeReqDispls, recvNodeReqDispls;
    int totalSendNodeReq = 0, totalRecvNodeReq = 0;
    prefixSums(sendNodeReqCounts, sendNodeReqDispls, totalSendNodeReq);
    prefixSums(recvNodeReqCounts, recvNodeReqDispls, totalRecvNodeReq);

    std::vector<int> sendNodeReqBuf(totalSendNodeReq);
    for (int r = 0; r < size; ++r) {
        int offset = sendNodeReqDispls[r];
        for (std::size_t k = 0; k < needNodes[r].size(); ++k) {
            sendNodeReqBuf[offset + static_cast<int>(k)] = needNodes[r][k];
        }
    }

    std::vector<int> recvNodeReqBuf(totalRecvNodeReq);
    MPI_Alltoallv(sendNodeReqBuf.data(),
                  sendNodeReqCounts.data(), sendNodeReqDispls.data(), MPI_INT,
                  recvNodeReqBuf.data(),
                  recvNodeReqCounts.data(), recvNodeReqDispls.data(), MPI_INT,
                  comm);

    // 4b. Phase 2: owners respond with coordinates in the same order

    std::vector<int> sendCoordCounts(size, 0);
    for (int r = 0; r < size; ++r) {
        sendCoordCounts[r] = recvNodeReqCounts[r] * nd;
    }

    std::vector<int> recvCoordCounts(size, 0);
    MPI_Alltoall(sendCoordCounts.data(), 1, MPI_INT,
                 recvCoordCounts.data(), 1, MPI_INT,
                 comm);

    std::vector<int> sendCoordDispls, recvCoordDispls;
    int totalSendCoord = 0, totalRecvCoord = 0;
    prefixSums(sendCoordCounts, sendCoordDispls, totalSendCoord);
    prefixSums(recvCoordCounts, recvCoordDispls, totalRecvCoord);

    std::vector<double> sendCoordBuf(totalSendCoord);

    // Fill sendCoordBuf: for each destination, in the order we received requests
    for (int src = 0; src < size; ++src) {
        int countNodesFromSrc = recvNodeReqCounts[src]; // nodes requested by 'src' from this rank
        int reqOffset         = recvNodeReqDispls[src];
        int coordOffset       = sendCoordDispls[src];

        for (int k = 0; k < countNodesFromSrc; ++k) {
            int gnode = recvNodeReqBuf[reqOffset + k];

            int ownerLocalStart = nodedist[rank];
            int ownerLocalEnd   = nodedist[rank + 1];

            if (gnode < ownerLocalStart || gnode >= ownerLocalEnd) {
                error("migrateMeshWithParMETIS: node request to wrong owner");
            }

            int lidx = gnode - ownerLocalStart; // local index in mesh_in.p

            for (int j = 0; j < nd; ++j) {
                sendCoordBuf[coordOffset + k * nd + j] =
                    mesh_in.p[j + lidx * nd];
            }
        }
    }

    std::vector<double> recvCoordBuf(totalRecvCoord);
    MPI_Alltoallv(sendCoordBuf.data(),
                  sendCoordCounts.data(), sendCoordDispls.data(), MPI_DOUBLE,
                  recvCoordBuf.data(),
                  recvCoordCounts.data(), recvCoordDispls.data(), MPI_DOUBLE,
                  comm);

    // ------------------------------------------------------------------
    // 5. BUILD FINAL LOCAL NODE LIST & GLOBAL→LOCAL MAPPING
    // ------------------------------------------------------------------

    mesh_out.nd  = nd;
    mesh_out.nve = nve;
    mesh_out.ne  = ne_local_out;
    mesh_out.np  = np_local_out;

    mesh_out.t.resize(static_cast<std::size_t>(nve) * ne_local_out);
    mesh_out.p.resize(static_cast<std::size_t>(nd)  * np_local_out);
    mesh_out.elemGlobalID.resize(ne_local_out);  // NEW

    // 5a. Build owner-based offsets into local node space
    std::vector<int> ownerNodeOffsets(size);
    {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            ownerNodeOffsets[r] = offset;
            offset += static_cast<int>(needNodes[r].size());
        }
    }

    // 5b. Global→local index map
    std::vector<int> g2l(np_global, -1);
    for (int r = 0; r < size; ++r) {
        int base = ownerNodeOffsets[r];
        for (std::size_t k = 0; k < needNodes[r].size(); ++k) {
            int gnode = needNodes[r][k];
            int lnode = base + static_cast<int>(k);
            g2l[gnode] = lnode;
        }
    }

    // 5c. Place received coordinates into mesh_out.p
    for (int owner = 0; owner < size; ++owner) {
        int numNodesFromOwner = sendNodeReqCounts[owner]; // nodes we requested from 'owner'
        if (numNodesFromOwner == 0) continue;

        int coordOffset = recvCoordDispls[owner]; // in doubles
        int baseLocal   = ownerNodeOffsets[owner];

        for (int k = 0; k < numNodesFromOwner; ++k) {
            int localNode = baseLocal + k;
            for (int j = 0; j < nd; ++j) {
                mesh_out.p[j + localNode * nd] =
                    recvCoordBuf[coordOffset + k * nd + j];
            }
        }
    }

    // ------------------------------------------------------------------
    // 6. REMAP CONNECTIVITY FROM GLOBAL NODE IDS → LOCAL NODE IDS
    // ------------------------------------------------------------------

    // Copy new element connectivity (still global IDs) into mesh_out.t,
    // and copy global element IDs into mesh_out.elemGlobalID
    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            mesh_out.t[j + e * nve] = recvElemConn[j + e * nve];
        }
        mesh_out.elemGlobalID[e] = recvElemGlobalID[e]; // NEW
    }

    // Replace global node IDs with local node IDs
    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            int gnode = mesh_out.t[j + e * nve];
            int lnode = g2l[gnode];
            if (lnode < 0) {
                error("migrateMeshWithParMETIS: missing local mapping for node");
            }
            mesh_out.t[j + e * nve] = lnode;
        }
    }
}

// // Parallel VTK reader: rank 0 reads full file, then distributes slices
// void readParMeshFromVTKFile(const std::string& filename,
//                             Mesh& mesh,
//                             MPI_Comm comm)
// {
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);
// 
//     Mesh globalMesh;   // only used on rank 0
// 
//     int nd_global  = 0;
//     int np_global  = 0;
//     int nve_global = 0;
//     int ne_global  = 0;
// 
//     // 1. Rank 0 reads the full mesh from VTK
//     if (rank == 0) {
//         readMeshFromVTKFile(filename, globalMesh);
//         nd_global  = globalMesh.nd;
//         np_global  = globalMesh.np;
//         nve_global = globalMesh.nve;
//         ne_global  = globalMesh.ne;
//     }
// 
//     // Broadcast global metadata
//     int header[4] = { nd_global, np_global, nve_global, ne_global };
//     MPI_Bcast(header, 4, MPI_INT, 0, comm);
//     nd_global  = header[0];
//     np_global  = header[1];
//     nve_global = header[2];
//     ne_global  = header[3];
// 
//     // 2. Compute local np, ne and offsets for this rank
//     int np_local   = 0;
//     int nodeOffset = 0;
//     int ne_local   = 0;
//     int elemOffset = 0;
// 
//     computeLocalRange(np_global, size, rank, np_local, nodeOffset);
//     computeLocalRange(ne_global, size, rank, ne_local, elemOffset);
// 
//     // 3. Initialize local mesh metadata and allocate storage
//     mesh.nd  = nd_global;
//     mesh.nve = nve_global;
//     mesh.np  = np_local;
//     mesh.ne  = ne_local;
// 
//     mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
//     mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);
// 
//     if (np_local == 0 && ne_local == 0) {
//         // This rank owns no nodes/elements; nothing to do
//         return;
//     }
// 
//     // 4. Distribute nodal coordinates (mesh.p), column-major: p(j + i*nd)
//     if (rank == 0) {
//         // Rank 0: fill its local block directly
//         for (int i = 0; i < np_local; ++i) {
//             int gnode = nodeOffset + i;
//             for (int j = 0; j < nd_global; ++j) {
//                 mesh.p[j + i * nd_global] =
//                     globalMesh.p[j + gnode * nd_global];
//             }
//         }
// 
//         // Send other ranks their blocks
//         for (int r = 1; r < size; ++r) {
//             int np_r, nodeOff_r;
//             computeLocalRange(np_global, size, r, np_r, nodeOff_r);
//             if (np_r == 0) continue;
// 
//             std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
//             for (int i = 0; i < np_r; ++i) {
//                 int gnode = nodeOff_r + i;
//                 for (int j = 0; j < nd_global; ++j) {
//                     sendP[j + i * nd_global] =
//                         globalMesh.p[j + gnode * nd_global];
//                 }
//             }
// 
//             MPI_Send(sendP.data(),
//                      static_cast<int>(sendP.size()),
//                      MPI_DOUBLE, r, 300, comm);
//         }
//     } else {
//         // Non-root ranks receive their p block
//         if (np_local > 0) {
//             std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
//             MPI_Recv(recvP.data(),
//                      static_cast<int>(recvP.size()),
//                      MPI_DOUBLE, 0, 300, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < np_local; ++i)
//                 for (int j = 0; j < nd_global; ++j)
//                     mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
//         }
//     }
// 
//     // 5. Distribute connectivity (mesh.t), column-major: t(j + i*nve)
//     if (rank == 0) {
//         // Rank 0: fill its local block
//         for (int i = 0; i < ne_local; ++i) {
//             int gele = elemOffset + i;
//             for (int j = 0; j < nve_global; ++j) {
//                 mesh.t[j + i * nve_global] =
//                     globalMesh.t[j + gele * nve_global];
//             }
//         }
// 
//         // Send other ranks their blocks
//         for (int r = 1; r < size; ++r) {
//             int ne_r, elemOff_r;
//             computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
//             if (ne_r == 0) continue;
// 
//             std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
//             for (int i = 0; i < ne_r; ++i) {
//                 int gele = elemOff_r + i;
//                 for (int j = 0; j < nve_global; ++j) {
//                     sendT[j + i * nve_global] =
//                         globalMesh.t[j + gele * nve_global];
//                 }
//             }
// 
//             MPI_Send(sendT.data(),
//                      static_cast<int>(sendT.size()),
//                      MPI_INT, r, 301, comm);
//         }
//     } else {
//         // Non-root ranks receive their t block
//         if (ne_local > 0) {
//             std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
//             MPI_Recv(recvT.data(),
//                      static_cast<int>(recvT.size()),
//                      MPI_INT, 0, 301, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < ne_local; ++i)
//                 for (int j = 0; j < nve_global; ++j)
//                     mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
//         }
//     }
// }
// 
// // Parallel VTU reader: rank 0 reads full VTU and distributes slices
// void readParMeshFromVTUFile(const std::string& filename,
//                             Mesh& mesh,
//                             MPI_Comm comm)
// {
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);
// 
//     Mesh globalMesh;   // only used on rank 0
// 
//     int nd_global  = 0;
//     int np_global  = 0;
//     int nve_global = 0;
//     int ne_global  = 0;
// 
//     // 1. Rank 0 reads full .vtu file using the serial reader
//     if (rank == 0) {
//         readMeshFromVTUFile(filename, globalMesh);
//         nd_global  = globalMesh.nd;
//         np_global  = globalMesh.np;
//         nve_global = globalMesh.nve;
//         ne_global  = globalMesh.ne;
//     }
// 
//     // Broadcast global metadata
//     int header[4] = { nd_global, np_global, nve_global, ne_global };
//     MPI_Bcast(header, 4, MPI_INT, 0, comm);
//     nd_global  = header[0];
//     np_global  = header[1];
//     nve_global = header[2];
//     ne_global  = header[3];
// 
//     // 2. Compute local np, ne and offsets for this rank
//     int np_local   = 0;
//     int nodeOffset = 0;
//     int ne_local   = 0;
//     int elemOffset = 0;
// 
//     computeLocalRange(np_global, size, rank, np_local, nodeOffset);
//     computeLocalRange(ne_global, size, rank, ne_local, elemOffset);
// 
//     // 3. Initialize local mesh metadata and allocate storage
//     mesh.nd  = nd_global;
//     mesh.nve = nve_global;
//     mesh.np  = np_local;
//     mesh.ne  = ne_local;
// 
//     mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
//     mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);
// 
//     if (np_local == 0 && ne_local == 0) {
//         // This rank owns nothing; done
//         return;
//     }
// 
//     // 4. Distribute nodal coordinates (mesh.p), column-major: p(j + i*nd)
//     if (rank == 0) {
//         // Rank 0: fill its local block directly
//         for (int i = 0; i < np_local; ++i) {
//             int gnode = nodeOffset + i;  // global node index
//             for (int j = 0; j < nd_global; ++j) {
//                 mesh.p[j + i * nd_global] =
//                     globalMesh.p[j + gnode * nd_global];
//             }
//         }
// 
//         // Send other ranks their blocks
//         for (int r = 1; r < size; ++r) {
//             int np_r, nodeOff_r;
//             computeLocalRange(np_global, size, r, np_r, nodeOff_r);
//             if (np_r == 0) continue;
// 
//             std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
//             for (int i = 0; i < np_r; ++i) {
//                 int gnode = nodeOff_r + i;
//                 for (int j = 0; j < nd_global; ++j) {
//                     sendP[j + i * nd_global] =
//                         globalMesh.p[j + gnode * nd_global];
//                 }
//             }
// 
//             MPI_Send(sendP.data(),
//                      static_cast<int>(sendP.size()),
//                      MPI_DOUBLE, r, 400, comm);
//         }
//     } else {
//         // Non-root ranks receive their p block
//         if (np_local > 0) {
//             std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
//             MPI_Recv(recvP.data(),
//                      static_cast<int>(recvP.size()),
//                      MPI_DOUBLE, 0, 400, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < np_local; ++i)
//                 for (int j = 0; j < nd_global; ++j)
//                     mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
//         }
//     }
// 
//     // 5. Distribute connectivity (mesh.t), column-major: t(j + i*nve)
//     if (rank == 0) {
//         // Rank 0: fill its local block
//         for (int i = 0; i < ne_local; ++i) {
//             int gele = elemOffset + i;   // global element index
//             for (int j = 0; j < nve_global; ++j) {
//                 mesh.t[j + i * nve_global] =
//                     globalMesh.t[j + gele * nve_global];
//             }
//         }
// 
//         // Send other ranks their blocks
//         for (int r = 1; r < size; ++r) {
//             int ne_r, elemOff_r;
//             computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
//             if (ne_r == 0) continue;
// 
//             std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
//             for (int i = 0; i < ne_r; ++i) {
//                 int gele = elemOff_r + i;
//                 for (int j = 0; j < nve_global; ++j) {
//                     sendT[j + i * nve_global] =
//                         globalMesh.t[j + gele * nve_global];
//                 }
//             }
// 
//             MPI_Send(sendT.data(),
//                      static_cast<int>(sendT.size()),
//                      MPI_INT, r, 401, comm);
//         }
//     } else {
//         // Non-root ranks receive their t block
//         if (ne_local > 0) {
//             std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
//             MPI_Recv(recvT.data(),
//                      static_cast<int>(recvT.size()),
//                      MPI_INT, 0, 401, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < ne_local; ++i)
//                 for (int j = 0; j < nve_global; ++j)
//                     mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
//         }
//     }
// }
// 
// // Parallel reader: rank 0 reads the .msh file, then distributes slices of p and t.
// void readParMeshFromMshV2File(const std::string& filename, Mesh& mesh, MPI_Comm comm)
// {
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);
// 
//     Mesh globalMesh;    // only used on rank 0
// 
//     int nd_global  = 0;
//     int np_global  = 0;
//     int nve_global = 0;
//     int ne_global  = 0;
// 
//     // ------------------------------------------------------------
//     // 1. Rank 0 reads the full mesh
//     // ------------------------------------------------------------
//     if (rank == 0) {
//         readMeshFromMshV2File(filename, globalMesh);
//         nd_global  = globalMesh.nd;
//         np_global  = globalMesh.np;
//         nve_global = globalMesh.nve;
//         ne_global  = globalMesh.ne;
//     }
// 
//     // Broadcast global sizes to all ranks
//     int header[4] = { nd_global, np_global, nve_global, ne_global };
//     MPI_Bcast(header, 4, MPI_INT, 0, comm);
//     nd_global  = header[0];
//     np_global  = header[1];
//     nve_global = header[2];
//     ne_global  = header[3];
// 
//     // ------------------------------------------------------------
//     // 2. Compute local np, ne and offsets for this rank
//     // ------------------------------------------------------------
//     int np_local   = 0;
//     int nodeOffset = 0;
//     int ne_local   = 0;
//     int elemOffset = 0;
// 
//     computeLocalRange(np_global, size, rank, np_local, nodeOffset);
//     computeLocalRange(ne_global, size, rank, ne_local, elemOffset);
// 
//     // ------------------------------------------------------------
//     // 3. Initialize local mesh metadata and allocate storage
//     // ------------------------------------------------------------
//     mesh.nd  = nd_global;   // same on all ranks
//     mesh.nve = nve_global;  // same on all ranks
//     mesh.np  = np_local;    // local
//     mesh.ne  = ne_local;    // local
// 
//     mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
//     mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);
// 
//     // Quick exit if this rank has no nodes and no elements
//     if (np_local == 0 && ne_local == 0) {
//         return;
//     }
// 
//     // ------------------------------------------------------------
//     // 4. Distribute nodes (mesh.p), column-major layout
//     //    p(j + i*nd), i = 0..np-1, j = 0..nd-1
//     // ------------------------------------------------------------
//     // Rank 0 sends each rank its contiguous block of nodes (by node index).
//     // Because of column-major layout, we pack into a temporary buffer per rank.
//     if (rank == 0) {
//         // For rank 0: fill its local mesh.p directly from globalMesh.p
//         for (int i = 0; i < np_local; ++i) {
//             int gnode = nodeOffset + i;  // global node index
//             for (int j = 0; j < nd_global; ++j) {
//                 mesh.p[j + i * nd_global] =
//                     globalMesh.p[j + gnode * nd_global];
//             }
//         }
// 
//         // For other ranks: pack and send
//         for (int r = 1; r < size; ++r) {
//             int np_r, nodeOff_r;
//             computeLocalRange(np_global, size, r, np_r, nodeOff_r);
// 
//             if (np_r == 0) continue;
// 
//             std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
//             for (int i = 0; i < np_r; ++i) {
//                 int gnode = nodeOff_r + i;
//                 for (int j = 0; j < nd_global; ++j) {
//                     sendP[j + i * nd_global] =
//                         globalMesh.p[j + gnode * nd_global];
//                 }
//             }
// 
//             MPI_Send(sendP.data(),
//                      static_cast<int>(sendP.size()),
//                      MPI_DOUBLE, r, 100, comm);
//         }
//     } else {
//         // Non-root ranks receive their p block
//         if (np_local > 0) {
//             std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
//             MPI_Recv(recvP.data(),
//                      static_cast<int>(recvP.size()),
//                      MPI_DOUBLE, 0, 100, comm, MPI_STATUS_IGNORE);
//             // Copy into mesh.p (same layout)
//             for (int i = 0; i < np_local; ++i)
//                 for (int j = 0; j < nd_global; ++j)
//                     mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
//         }
//     }
// 
//     // ------------------------------------------------------------
//     // 5. Distribute elements (mesh.t), column-major layout
//     //    t(j + i*nve), i = 0..ne-1, j = 0..nve-1
//     // ------------------------------------------------------------
//     if (rank == 0) {
//         // Rank 0: fill its local mesh.t
//         for (int i = 0; i < ne_local; ++i) {
//             int gele = elemOffset + i;  // global element index
//             for (int j = 0; j < nve_global; ++j) {
//                 mesh.t[j + i * nve_global] =
//                     globalMesh.t[j + gele * nve_global];
//             }
//         }
// 
//         // Other ranks: pack and send
//         for (int r = 1; r < size; ++r) {
//             int ne_r, elemOff_r;
//             computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
//             if (ne_r == 0) continue;
// 
//             std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
//             for (int i = 0; i < ne_r; ++i) {
//                 int gele = elemOff_r + i;
//                 for (int j = 0; j < nve_global; ++j) {
//                     sendT[j + i * nve_global] =
//                         globalMesh.t[j + gele * nve_global];
//                 }
//             }
// 
//             MPI_Send(sendT.data(),
//                      static_cast<int>(sendT.size()),
//                      MPI_INT, r, 101, comm);
//         }
//     } else {
//         // Non-root ranks receive their t block
//         if (ne_local > 0) {
//             std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
//             MPI_Recv(recvT.data(),
//                      static_cast<int>(recvT.size()),
//                      MPI_INT, 0, 101, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < ne_local; ++i)
//                 for (int j = 0; j < nve_global; ++j)
//                     mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
//         }
//     }
// }
// 
// // Parallel reader: rank 0 reads gmsh v4 file, distributes slices to all ranks
// void readParMeshFromMshV4File(const std::string& filename,
//                               Mesh& mesh,
//                               MPI_Comm comm)
// {
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);
// 
//     Mesh globalMesh;   // only used on rank 0
// 
//     int nd_global  = 0;
//     int np_global  = 0;
//     int nve_global = 0;
//     int ne_global  = 0;
// 
//     // 1. Rank 0 reads the full mesh
//     if (rank == 0) {
//         readMeshFromMshV4File(filename, globalMesh);
//         nd_global  = globalMesh.nd;
//         np_global  = globalMesh.np;
//         nve_global = globalMesh.nve;
//         ne_global  = globalMesh.ne;
//     }
// 
//     // Broadcast global metadata
//     int header[4] = { nd_global, np_global, nve_global, ne_global };
//     MPI_Bcast(header, 4, MPI_INT, 0, comm);
//     nd_global  = header[0];
//     np_global  = header[1];
//     nve_global = header[2];
//     ne_global  = header[3];
// 
//     // 2. Compute local np, ne and offsets
//     int np_local   = 0;
//     int nodeOffset = 0;
//     int ne_local   = 0;
//     int elemOffset = 0;
// 
//     computeLocalRange(np_global, size, rank, np_local, nodeOffset);
//     computeLocalRange(ne_global, size, rank, ne_local, elemOffset);
// 
//     // 3. Initialize local mesh metadata and allocate
//     mesh.nd  = nd_global;
//     mesh.nve = nve_global;
//     mesh.np  = np_local;
//     mesh.ne  = ne_local;
// 
//     mesh.p.resize(static_cast<std::size_t>(mesh.nd) * mesh.np);
//     mesh.t.resize(static_cast<std::size_t>(mesh.nve) * mesh.ne);
// 
//     if (np_local == 0 && ne_local == 0) {
//         // This rank owns nothing; done
//         return;
//     }
// 
//     // 4. Distribute nodes (mesh.p), column-major: p(j + i*nd)
//     if (rank == 0) {
//         // Fill rank 0's local p directly
//         for (int i = 0; i < np_local; ++i) {
//             int gnode = nodeOffset + i;
//             for (int j = 0; j < nd_global; ++j) {
//                 mesh.p[j + i * nd_global] =
//                     globalMesh.p[j + gnode * nd_global];
//             }
//         }
// 
//         // Send to other ranks
//         for (int r = 1; r < size; ++r) {
//             int np_r, nodeOff_r;
//             computeLocalRange(np_global, size, r, np_r, nodeOff_r);
//             if (np_r == 0) continue;
// 
//             std::vector<double> sendP(static_cast<std::size_t>(np_r) * nd_global);
//             for (int i = 0; i < np_r; ++i) {
//                 int gnode = nodeOff_r + i;
//                 for (int j = 0; j < nd_global; ++j) {
//                     sendP[j + i * nd_global] =
//                         globalMesh.p[j + gnode * nd_global];
//                 }
//             }
// 
//             MPI_Send(sendP.data(),
//                      static_cast<int>(sendP.size()),
//                      MPI_DOUBLE, r, 200, comm);
//         }
//     } else {
//         // Receive p on non-root ranks
//         if (np_local > 0) {
//             std::vector<double> recvP(static_cast<std::size_t>(np_local) * nd_global);
//             MPI_Recv(recvP.data(),
//                      static_cast<int>(recvP.size()),
//                      MPI_DOUBLE, 0, 200, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < np_local; ++i)
//                 for (int j = 0; j < nd_global; ++j)
//                     mesh.p[j + i * nd_global] = recvP[j + i * nd_global];
//         }
//     }
// 
//     // 5. Distribute elements (mesh.t), column-major: t(j + i*nve)
//     if (rank == 0) {
//         // Rank 0 local elements
//         for (int i = 0; i < ne_local; ++i) {
//             int gele = elemOffset + i;
//             for (int j = 0; j < nve_global; ++j) {
//                 mesh.t[j + i * nve_global] =
//                     globalMesh.t[j + gele * nve_global];
//             }
//         }
// 
//         // Send to other ranks
//         for (int r = 1; r < size; ++r) {
//             int ne_r, elemOff_r;
//             computeLocalRange(ne_global, size, r, ne_r, elemOff_r);
//             if (ne_r == 0) continue;
// 
//             std::vector<int> sendT(static_cast<std::size_t>(ne_r) * nve_global);
//             for (int i = 0; i < ne_r; ++i) {
//                 int gele = elemOff_r + i;
//                 for (int j = 0; j < nve_global; ++j) {
//                     sendT[j + i * nve_global] =
//                         globalMesh.t[j + gele * nve_global];
//                 }
//             }
// 
//             MPI_Send(sendT.data(),
//                      static_cast<int>(sendT.size()),
//                      MPI_INT, r, 201, comm);
//         }
//     } else {
//         // Receive t on non-root ranks
//         if (ne_local > 0) {
//             std::vector<int> recvT(static_cast<std::size_t>(ne_local) * nve_global);
//             MPI_Recv(recvT.data(),
//                      static_cast<int>(recvT.size()),
//                      MPI_INT, 0, 201, comm, MPI_STATUS_IGNORE);
// 
//             for (int i = 0; i < ne_local; ++i)
//                 for (int j = 0; j < nve_global; ++j)
//                     mesh.t[j + i * nve_global] = recvT[j + i * nve_global];
//         }
//     }
// }
