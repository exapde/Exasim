#include <mpi.h>
#include <parmetis.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>

/// Parallel mesh partitioning using ParMETIS_V3_PartMeshKway
///
/// elmdist: global element distribution over MPI ranks, size = nprocs + 1.
/// eind:    local element-to-node connectivity (flattened, global node IDs)
///
/// Output:
///   epart_local[e] = partition ID of local element e
///
void partitionMeshParMETIS(std::vector<idx_t>       &epart_local,
                           std::vector<idx_t>       &eind,      // local connectivity
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

    if (ne_local > 0 && eind.size() != (size_t)ne_local * nve) {
        if (rank == 0) {
            std::cerr << "partitionMeshParMETIS: eind.size() = " << eind.size()
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
    idx_t *eind_ptr    = eind.data();

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
