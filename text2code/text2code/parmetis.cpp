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

