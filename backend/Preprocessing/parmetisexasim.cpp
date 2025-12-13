#ifndef __PARMETISEXASIM
#define __PARMETISEXASIM

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
        std::cout << "Finished partitioning mesh using ParMETIS (edgecut = " << edgecut << ")" << std::endl;
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
    (void)ne_global; // unused but kept for clarity

    const int gstart = elmdist[rank];
    std::vector<int> elemGlobalID_in(ne_local_in);
    for (int e = 0; e < ne_local_in; ++e) {
        elemGlobalID_in[e] = gstart + e;
    }

    // ------------------------------------------------------------------
    // 1. ELEMENT MIGRATION (based on ParMETIS partition epart_local)
    // ------------------------------------------------------------------
    std::vector<int> sendElemCounts(size, 0);
    for (int e = 0; e < ne_local_in; ++e) {
        int dest = static_cast<int>(epart_local[e]);
        if (dest < 0 || dest >= size) {
            error("migrateMeshWithParMETIS: epart_local has invalid destination rank");
        }
        sendElemCounts[dest] += 1;
    }

    std::vector<int> recvElemCounts(size, 0);
    MPI_Alltoall(sendElemCounts.data(), 1, MPI_INT,
                 recvElemCounts.data(), 1, MPI_INT,
                 comm);

    std::vector<int> sendElemDispls, recvElemDispls;
    int totalSendElems = 0, totalRecvElems = 0;
    prefixSums(sendElemCounts, sendElemDispls, totalSendElems);
    prefixSums(recvElemCounts, recvElemDispls, totalRecvElems);

    std::vector<int> sendElemConn(static_cast<std::size_t>(totalSendElems) * nve);
    std::vector<int> sendElemGlobalID(static_cast<std::size_t>(totalSendElems));

    std::vector<int> elemOffsetPerDest(size, 0);

    for (int e = 0; e < ne_local_in; ++e) {
        int dest           = static_cast<int>(epart_local[e]);
        int localIdxInDest = elemOffsetPerDest[dest]++;
        int elemPos        = sendElemDispls[dest] + localIdxInDest;

        for (int j = 0; j < nve; ++j) {
            int gnode = mesh_in.t[j + e * nve]; // currently global node ID
            sendElemConn[j + elemPos * nve] = gnode;
        }

        sendElemGlobalID[elemPos] = elemGlobalID_in[e];
    }

    
    std::vector<int> recvElemConn(static_cast<std::size_t>(totalRecvElems) * nve);
    std::vector<int> recvElemGlobalID(static_cast<std::size_t>(totalRecvElems));

    // 1e. Prepare counts/displs in *ints* for connectivity exchange
    std::vector<int> sendElemCountsInts(size), recvElemCountsInts(size);
    std::vector<int> sendElemDisplsInts(size), recvElemDisplsInts(size);
    
    for (int r = 0; r < size; ++r) {
        sendElemCountsInts[r] = sendElemCounts[r] * nve;
        recvElemCountsInts[r] = recvElemCounts[r] * nve;
        sendElemDisplsInts[r] = sendElemDispls[r] * nve;
        recvElemDisplsInts[r] = recvElemDispls[r] * nve;
    }


    MPI_Alltoallv(sendElemConn.data(),
                  sendElemCountsInts.data(), sendElemDisplsInts.data(), MPI_INT,
                  recvElemConn.data(),
                  recvElemCountsInts.data(), recvElemDisplsInts.data(), MPI_INT,
                  comm);

    MPI_Alltoallv(sendElemGlobalID.data(),
                  sendElemCounts.data(),
                  sendElemDispls.data(),   // in "elements"
                  MPI_INT,
                  recvElemGlobalID.data(),
                  recvElemCounts.data(),
                  recvElemDispls.data(),
                  MPI_INT,
                  comm);

    const int ne_local_out = totalRecvElems;

    // ------------------------------------------------------------------
    // 2. BUILD NODE OWNERSHIP INFORMATION (contiguous distribution)
    // ------------------------------------------------------------------
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

    auto findNodeOwner = [&](int gnode) -> int {
        for (int r = 0; r < size; ++r) {
            if (gnode >= nodedist[r] && gnode < nodedist[r + 1]) {
                return r;
            }
        }
        return -1;
    };

    // ------------------------------------------------------------------
    // 3. DETERMINE WHICH GLOBAL NODES THIS RANK NEEDS
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

    std::vector<std::vector<int>> needNodes(size);
    for (int g = 0; g < np_global; ++g) {
        if (nodeMask[g]) {
            int owner = findNodeOwner(g);
            if (owner < 0) {
                error("migrateMeshWithParMETIS: could not find owner for global node");
            }
            needNodes[owner].push_back(g);
        }
    }

    int np_local_out = 0;
    for (int r = 0; r < size; ++r) {
        np_local_out += static_cast<int>(needNodes[r].size());
    }

    // ------------------------------------------------------------------
    // 4. REQUEST NODE COORDINATES FROM OWNERS
    // ------------------------------------------------------------------
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

    for (int src = 0; src < size; ++src) {
        int countNodesFromSrc = recvNodeReqCounts[src];
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
    mesh_out.elemGlobalID.resize(ne_local_out);
    mesh_out.nodeGlobalID.resize(np_local_out);   // NEW

    // 5a. Owner-based offsets into local node space
    std::vector<int> ownerNodeOffsets(size);
    {
        int offset = 0;
        for (int r = 0; r < size; ++r) {
            ownerNodeOffsets[r] = offset;
            offset += static_cast<int>(needNodes[r].size());
        }
    }

    // 5b. Global→local map
    std::vector<int> g2l(np_global, -1);
    for (int r = 0; r < size; ++r) {
        int base = ownerNodeOffsets[r];
        for (std::size_t k = 0; k < needNodes[r].size(); ++k) {
            int gnode = needNodes[r][k];
            int lnode = base + static_cast<int>(k);
            g2l[gnode] = lnode;

            // Fill nodeGlobalID here: local node -> its global ID
            mesh_out.nodeGlobalID[lnode] = gnode;     // NEW
        }
    }

    // 5c. Place received coordinates into mesh_out.p
    for (int owner = 0; owner < size; ++owner) {
        int numNodesFromOwner = sendNodeReqCounts[owner];
        if (numNodesFromOwner == 0) continue;

        int coordOffset = recvCoordDispls[owner];
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
    for (int e = 0; e < ne_local_out; ++e) {
        for (int j = 0; j < nve; ++j) {
            mesh_out.t[j + e * nve] = recvElemConn[j + e * nve];
        }
        mesh_out.elemGlobalID[e] = recvElemGlobalID[e];
    }

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

    if (rank == 0) {
        std::cout << "Finished migrating Mesh after ParMETIS (ne_local_out = " << ne_local_out << ")" << std::endl;
    }  
}

// // ---- Build e2e: for each element 'e' and local face 'lf', write neighbor element id or -1 ----
// // e2e must be sized to ne * nfe (ints).
// int mke2e_hash(int* e2e,
//                const int* e2n,          // [ne * nne] global node ids per element
//                const int* local_faces,  // [nfe * nnf] local node indices per local face
//                int ne, int nne, int nnf, int nfe)
// {
//     using Val = std::pair<int,int>; // (elem, lf) of the first owner of the face
// 
//     // init all neighbors to -1 (assume boundary until proven otherwise)
//     for (int i = 0; i < ne * nfe; ++i) e2e[i] = -1;
// 
//     // map holds faces seen once (awaiting their twin)
//     std::unordered_map<FaceKey, Val, FaceKeyHash> map;
//     map.reserve(static_cast<size_t>(ne) * static_cast<size_t>(nfe) * 13 / 10); // slack to keep LF low
// 
//     int out = 0;
// 
//     auto make_key = [&](int e, int lf) {
//         FaceKey k; k.len = static_cast<uint8_t>(nnf);
//         // gather global node ids of this face
//         for (int i = 0; i < nnf; ++i) {
//             const int ln = local_faces[lf * nnf + i];   // local node id on this face
//             k.a[i] = e2n[e * nne + ln];                 // map to global node id
//         }
//         // canonicalize (orientation-independent)
//         std::sort(k.a.begin(), k.a.begin() + nnf);
//         return k;
//     };
// 
//     for (int e = 0; e < ne; ++e) {
//         for (int lf = 0; lf < nfe; ++lf) {
//             FaceKey k = make_key(e, lf);
//             auto it = map.find(k);
// 
//             if (it == map.end()) {
//                 // first time seeing this face
//                 map.emplace(std::move(k), Val{e, lf});
//             } else {
//                 // found twin: set both directions and remove from map
//                 const auto [e0, lf0] = it->second;
// 
//                 e2e[e0 * nfe + lf0] = e;   // neighbor of the first owner is current element
//                 e2e[e  * nfe + lf ] = e0;  // neighbor of current element is the first owner
// 
//                 map.erase(it);
//                 ++out;
//             }
//         }
//     }
// 
//     return out; 
// }

void mke2e_global(int* e2e,                 // [ne * nfe], output GLOBAL neighbor IDs
                    const int* e2n,         // [ne * nne], local node IDs
                    const int* local_faces, // [nfe * nnf]
                    const int* elemGlobalID,// [ne]        global element IDs
                    int ne, int nne, int nnf, int nfe, int rank)
{
    
    // printf("%d %d %d %d %d\n", rank, ne, nne, nnf, nfe);
    // if (rank==0) {
    //   print2iarray(local_faces, nnf, nfe);
    //   print2iarray(e2n, nne, ne);
    //   print2iarray(elemGlobalID, 1, ne);
    // }
    // MPI_Barrier(EXASIM_COMM_WORLD);
    // error("here");

    // 1. Build local neighbor connectivity
    std::vector<int> e2e_local(ne * nfe);
    mke2e_hash(e2e_local.data(), e2n, local_faces, ne, nne, nnf, nfe);

    // 2. Convert local neighbor index → global neighbor ID
    for (int e = 0; e < ne; ++e) {
        for (int lf = 0; lf < nfe; ++lf) {
            int locNbr = e2e_local[e * nfe + lf];
            if (locNbr < 0) {
                e2e[e * nfe + lf] = -1;         // no neighbor / off-rank
            } else {
                e2e[e * nfe + lf] = elemGlobalID[locNbr];
            }
        }
    }

    if (rank == 0) {
        std::cout << "Finished making element-to-element connectivities on each subdomain" << std::endl;
    }    
}

void callParMetis(Mesh& mesh, const PDE& pde, MPI_Comm comm) 
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  
    mesh.elmdist = buildElmdistFromLocalCount(mesh.ne, comm);      
    
    if (pde.partitionfile == "")  {
      partitionMeshParMETIS(mesh.epart_local, mesh.t, mesh.elmdist, mesh.nve, mesh.nvf, size, comm);
    }
    else {        
      vector<double> tm;
      readarrayfrombinaryfile(make_path(pde.datapath, pde.partitionfile), tm);      
      int n1 = mesh.elmdist[rank];
      int n2 = mesh.elmdist[rank+1];
      int n = n2 - n1;
      mesh.epart_local.resize(n);
      for (int j = 0; j < n; j++)
        mesh.epart_local[j] = (idx_t) (tm[n1+j]-1);   
    }

    Mesh mesh_in;
    mesh_in.nd = mesh.nd;
    mesh_in.nve = mesh.nve;
    mesh_in.ne = mesh.ne;
    mesh_in.np = mesh.np;
    mesh_in.t.resize(mesh.nve*mesh.ne);
    mesh_in.p.resize(mesh.nd*mesh.np);
    for (int i = 0; i < mesh.nve*mesh.ne; i++) mesh_in.t[i] = mesh.t[i];
    for (int i = 0; i < mesh.nd*mesh.np; i++) mesh_in.p[i] = mesh.p[i];

    migrateMeshWithParMETIS(mesh_in, mesh.epart_local, mesh, comm);      

    mesh.localfaces.resize(mesh.nvf * mesh.nfe);           
    getelemface(mesh.localfaces.data(), mesh.dim, mesh.elemtype);    

    //std::vector<int> t2t_local(mesh.nfe*mesh.ne);      
    mesh.t2t.resize(mesh.nfe*mesh.ne);
    mke2e_global(mesh.t2t.data(), mesh.t.data(), mesh.localfaces.data(), 
                 mesh.elemGlobalID.data(), mesh.ne, mesh.nve, mesh.nvf, mesh.nfe, rank);     
}

// After mke2e_global: fill remote neighbors (across MPI ranks) in e2e.
//
// e2e          : [ne * nfe], input/output, GLOBAL neighbor element IDs
//                - already has local neighbors;
//                - entries == -1 are either boundary or cross-rank faces.
// e2n          : [ne * nne], local node IDs per element
// local_faces  : [nfe * nnf], local node indices per local face
// elemGlobalID : [ne], global element ID for each local element
// nodeGlobalID : [np], global node ID for each local node
// ne, nne      : #elements, #nodes per element
// nnf, nfe     : #nodes per face, #faces per element
// nbinfo: output array of 6*nbinfoSize ints with records:
//   (eLoc, lf, self_rank, eglonb, lfnb, neighbor_rank)
void mke2e_fill_first_neighbors(int*       e2e,
                                 const int* e2n,
                                 const int* local_faces,
                                 const int* elemGlobalID,
                                 const int* nodeGlobalID,
                                 int        ne,
                                 int        nne,
                                 int        nnf,
                                 int        nfe,
                                 MPI_Comm   comm,
                                 vector<int>  &nbinfoVec)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (nnf > 4) {
        error("mke2e_fill_remote_neighbors: nnf > 4 not supported in this implementation");
    }

    // ------------------------------------------------------------------
    // 1. Collect all faces with e2e == -1 (potential boundary or remote)
    //    Owner rule: owner = min(global_node_ids) % size
    // ------------------------------------------------------------------
    struct LocalFaceRecord {
        std::array<int,4> gnodes; // global node IDs of this face
        int               elemGlobal;
        int               elemLocal;
        int               lf;
    };

    std::vector<int>             sendFaceCounts(size, 0);
    std::vector<LocalFaceRecord> exportFaces;
    exportFaces.reserve(ne * nfe);

    for (int e = 0; e < ne; ++e) {
        int eGlob = elemGlobalID[e];
        for (int lf = 0; lf < nfe; ++lf) {
            if (e2e[e * nfe + lf] != -1) {
                // already has on-rank neighbor, skip
                continue;
            }

            LocalFaceRecord rec;
            rec.elemGlobal = eGlob;
            rec.elemLocal  = e;
            rec.lf         = lf;

            int minNode = std::numeric_limits<int>::max();
            for (int i = 0; i < nnf; ++i) {
                int ln  = local_faces[lf * nnf + i];  // local node in element
                int ln2 = e2n[e * nne + ln];          // local node in mesh
                int gn  = nodeGlobalID[ln2];          // global node ID
                rec.gnodes[i] = gn;
                if (gn < minNode) minNode = gn;
            }

            int owner = (minNode % size + size) % size;
            sendFaceCounts[owner] += 1;
            exportFaces.push_back(rec);
        }
    }

    // ------------------------------------------------------------------
    // 2. Pack faces into a send buffer and Alltoallv to owners
    //
    // Face record layout (ints):
    //   [0..nnf-1] : global node IDs
    //   [nnf]      : elemGlobal
    //   [nnf+1]    : elemLocal
    //   [nnf+2]    : lf
    // recordSize = nnf + 3
    // ------------------------------------------------------------------
    const int recordSize = nnf + 3;

    std::vector<int> sendCountsInts(size), sendDisplsInts;
    std::vector<int> recvCountsInts(size), recvDisplsInts;
    for (int r = 0; r < size; ++r) {
        sendCountsInts[r] = sendFaceCounts[r] * recordSize;
    }

    int totalSendInts = 0, totalRecvInts = 0;
    prefixSums(sendCountsInts, sendDisplsInts, totalSendInts);

    MPI_Alltoall(sendCountsInts.data(), 1, MPI_INT,
                 recvCountsInts.data(), 1, MPI_INT,
                 comm);

    prefixSums(recvCountsInts, recvDisplsInts, totalRecvInts);

    std::vector<int> sendBuf(totalSendInts);
    std::vector<int> faceOffsetPerDest(size, 0); // in faces

    for (const auto& rec : exportFaces) {
        int minNode = std::numeric_limits<int>::max();
        for (int i = 0; i < nnf; ++i)
            if (rec.gnodes[i] < minNode) minNode = rec.gnodes[i];

        int owner = (minNode % size + size) % size;

        int faceIdxInDest = faceOffsetPerDest[owner]++;
        int base = sendDisplsInts[owner] + faceIdxInDest * recordSize;

        for (int i = 0; i < nnf; ++i) {
            sendBuf[base + i] = rec.gnodes[i];
        }
        sendBuf[base + nnf    ] = rec.elemGlobal;
        sendBuf[base + nnf + 1] = rec.elemLocal;
        sendBuf[base + nnf + 2] = rec.lf;
    }

    std::vector<int> recvBuf(totalRecvInts);

    MPI_Alltoallv(sendBuf.data(),
                  sendCountsInts.data(), sendDisplsInts.data(), MPI_INT,
                  recvBuf.data(),
                  recvCountsInts.data(), recvDisplsInts.data(), MPI_INT,
                  comm);

    // ------------------------------------------------------------------
    // 3. On owner ranks: match faces by FaceKey and build neighbor info
    // ------------------------------------------------------------------
    struct FaceOwnerRec {
        int rank;        // original rank
        int elemGlobal;  // global element ID
        int elemLocal;   // local element on that rank
        int lf;          // local face index on that element
    };

    std::unordered_map<FaceKey, std::vector<FaceOwnerRec>, FaceKeyHash> faceMap;
    faceMap.reserve(static_cast<std::size_t>(totalRecvInts / recordSize) * 2);

    // Decode recvBuf into faceMap
    for (int src = 0; src < size; ++src) {
        int countInts = recvCountsInts[src];
        int baseInts  = recvDisplsInts[src];

        int numFaces = (recordSize > 0) ? (countInts / recordSize) : 0;
        for (int f = 0; f < numFaces; ++f) {
            int base = baseInts + f * recordSize;

            FaceOwnerRec rec;
            rec.rank       = src;
            rec.elemGlobal = recvBuf[base + nnf];
            rec.elemLocal  = recvBuf[base + nnf + 1];
            rec.lf         = recvBuf[base + nnf + 2];

            FaceKey key;
            key.len = static_cast<uint8_t>(nnf);
            for (int i = 0; i < nnf; ++i) {
                key.a[i] = recvBuf[base + i];
            }
            std::sort(key.a.begin(), key.a.begin() + nnf);

            auto& vec = faceMap[key];
            vec.push_back(rec);
        }
    }

    // Updates to send back: for each record, we want:
    //   selfElemLocal, selfLf, neighborElemGlobal, neighborLf, neighborRank
    struct Update {
        int selfElemLocal;
        int selfLf;
        int neighGlobal;
        int neighLf;
        int neighRank;
    };

    std::vector<std::vector<Update>> updates(size);

    for (auto& kv : faceMap) {
        auto& owners = kv.second;
        if (owners.size() == 2) {
            const auto& a = owners[0];
            const auto& b = owners[1];

            // To rank a.rank: self = a, neighbor = b
            updates[a.rank].push_back(
                Update{a.elemLocal, a.lf, b.elemGlobal, b.lf, b.rank});

            // To rank b.rank: self = b, neighbor = a
            updates[b.rank].push_back(
                Update{b.elemLocal, b.lf, a.elemGlobal, a.lf, a.rank});
        }
        // owners.size() == 1 → physical boundary
        // owners.size() > 2 → invalid non-conforming mesh
    }

    // ------------------------------------------------------------------
    // 4. Send updates back, fill e2e and nbinfo
    //
    // Each update record (ints):
    //   [0] = selfElemLocal
    //   [1] = selfLf
    //   [2] = neighGlobal
    //   [3] = neighLf
    //   [4] = neighRank
    //
    // updRecordSize = 5
    // ------------------------------------------------------------------
    const int updRecordSize = 5;

    std::vector<int> updSendCountsInts(size), updSendDisplsInts;
    std::vector<int> updRecvCountsInts(size), updRecvDisplsInts;

    for (int r = 0; r < size; ++r) {
        updSendCountsInts[r] = static_cast<int>(updates[r].size()) * updRecordSize;
    }

    int totalUpdSendInts = 0, totalUpdRecvInts = 0;
    prefixSums(updSendCountsInts, updSendDisplsInts, totalUpdSendInts);

    MPI_Alltoall(updSendCountsInts.data(), 1, MPI_INT,
                 updRecvCountsInts.data(), 1, MPI_INT,
                 comm);

    prefixSums(updRecvCountsInts, updRecvDisplsInts, totalUpdRecvInts);

    std::vector<int> updSendBuf(totalUpdSendInts);
    std::vector<int> updRecvBuf(totalUpdRecvInts);

    // Pack updates
    for (int r = 0; r < size; ++r) {
        int base = updSendDisplsInts[r];
        const auto& ur = updates[r];
        for (std::size_t k = 0; k < ur.size(); ++k) {
            updSendBuf[base + k * updRecordSize + 0] = ur[k].selfElemLocal;
            updSendBuf[base + k * updRecordSize + 1] = ur[k].selfLf;
            updSendBuf[base + k * updRecordSize + 2] = ur[k].neighGlobal;
            updSendBuf[base + k * updRecordSize + 3] = ur[k].neighLf;
            updSendBuf[base + k * updRecordSize + 4] = ur[k].neighRank;
        }
    }

    MPI_Alltoallv(updSendBuf.data(),
                  updSendCountsInts.data(), updSendDisplsInts.data(), MPI_INT,
                  updRecvBuf.data(),
                  updRecvCountsInts.data(), updRecvDisplsInts.data(), MPI_INT,
                  comm);

    // Collect neighbor info triples (6 ints per record)
    //std::vector<int> nbinfoVec;

    // Apply updates locally and build nbinfo
    for (int src = 0; src < size; ++src) {
        int countInts = updRecvCountsInts[src];
        int baseInts  = updRecvDisplsInts[src];
        int numUpd    = (updRecordSize > 0) ? (countInts / updRecordSize) : 0;

        for (int k = 0; k < numUpd; ++k) {
            int base      = baseInts + k * updRecordSize;
            int eLoc      = updRecvBuf[base + 0];
            int lf        = updRecvBuf[base + 1];
            int nghGlobal = updRecvBuf[base + 2];
            int lfnb      = updRecvBuf[base + 3];
            int nghRank   = updRecvBuf[base + 4];

            // overwrite -1 with neighbor global element ID
            e2e[eLoc * nfe + lf] = nghGlobal;

            // record (eLoc, lf, self_rank, eglonb, lfnb, neighbor_rank)
            nbinfoVec.insert(nbinfoVec.end(), {eLoc, lf, rank, nghGlobal, lfnb, nghRank});
            // nbinfoVec.push_back(eLoc);
            // nbinfoVec.push_back(lf);
            // nbinfoVec.push_back(rank);     // self_rank
            // nbinfoVec.push_back(nghGlobal);
            // nbinfoVec.push_back(lfnb);
            // nbinfoVec.push_back(nghRank);
        }
    }

    if (rank == 0) {
        std::cout << "Finished making first-neighbor connectivities between subdomains" << std::endl;
    }    
}

// Label boundary faces in t2t:
// For any face (l,e) with t2t[l + nfe*e] == -1 and matching boundary
// expression index k, set t2t[l + nfe*e] = -(k+1).
//
// - t2t: [nfe x ne], element-to-element; boundary faces initially have -1
// - t  : [nve x ne], element-to-node connectivity
// - localfaces: [nvf x nfe], local face connectivity
// - p  : [dim x np], nodal coordinates
// - bndexpr: array of boundary expressions, e.g., {"x*x + y*y <= 1"}
// - nbndexpr: number of boundary expressions
//
// Assumes eval_expr(char* expr, const double* x, int dim) is available.
int setboundaryfaces(
    int* t2t,                  // [nfe x ne], element-to-element; boundary faces initially -1
    const int* e2n,            // [nve x ne], element-to-node connectivity
    const int* localfaces,     // [nvf x nfe], local face connectivity
    const double* p,           // [dim x np], nodal coordinates
    char** bndexpr,            // boundary expressions
    int dim, int nve, int nvf, int nfe, int ne, int nbndexpr)
{
    int ind2 = (nvf == 1) ? 0 : 1;
    int count_bfaces = 0;

    // Loop over all elements and local faces
    for (int e = 0; e < ne; ++e) {
        for (int l = 0; l < nfe; ++l) {
            int face_pos = l + nfe * e;

            // Only process boundary faces (originally marked -1)
            if (t2t[face_pos] != -1)
                continue;

            bool labeled = false;

            // for (int i = 0; i < nvf; i++) {              
            //   int local_node = localfaces[i + l * nvf];
            //   int node = e2n[local_node + nve * e];  // global node
            //   const double* x = &p[dim * node];
            //   for (int k = 0; k < nbndexpr; ++k)
            //     if (eval_expr(bndexpr[k], x, dim)) 
            //       printf("%d %d %d %d %g %g\n", e, l, node, k, x[0], x[1]);  
            // }

            // Try matching boundary expressions
            for (int k = 0; k < nbndexpr; ++k) {
                bool ok = true;
                int check_pts[3] = {0, ind2, nvf - 1};

                for (int m = 0; m < 3; ++m) {
                    int i = check_pts[m];

                    // Local node index on face l
                    int local_node = localfaces[i + l * nvf];
                    int node = e2n[local_node + nve * e];  // global node
                    const double* x = &p[dim * node];

                    if (!eval_expr(bndexpr[k], x, dim)) {
                        ok = false;
                        break;
                    }
                }

                if (ok) {
                    t2t[face_pos] = -(k + 1);  // encode boundary index
                    labeled = true;
                    break;
                }
            }

            if (labeled)
                count_bfaces++;
        }
    }

    return count_bfaces;
}

// t2t: [nfe x ne], element-to-element (global neighbor element IDs).
//      Boundary faces have been labeled by labelboundaryfaces_in_t2t as
//      t2t = -(k+1), where k is the boundary-expression index.
// t:   [nve x ne], element-to-node connectivity (local node indices into p)
// localfaces: [nvf x nfe], local node ids per local face (0..nve-1)
// p:   [dim x np], nodal coordinates for local nodes
// elemGlobalID: [ne], global element ID for each local element
// prd_f1, prd_f2: [nprd], periodic boundary pairs; each entry is 1-based
//                 index into boundary expressions (like setperiodicfaces).
// expr1, expr2: [nprd * ncomp], mapping expressions for each periodic pair.
// dim, nve, nvf, nfe, ne: mesh sizes
// nprd: number of periodic boundary pairs
// ncomp: number of mapping components (for xiny)
// comm: MPI communicator
// nbinfoVec: will be appended with records
//            (eLoc, lf, self_rank, eglonb, lfnb, neighbor_rank)
//            for cross-rank periodic neighbors.
void setperiodicfaces(
    int*       t2t,
    const int* t,
    const int* localfaces,
    const double* p,
    const int* elemGlobalID,
    const int* prd_f1,
    const int* prd_f2,
    char** expr1,
    char** expr2,
    int dim,
    int nve,
    int nvf,
    int nfe,
    int ne,
    int nprd,
    int ncomp,
    int nbfaces,
    MPI_Comm comm,
    std::vector<int>& nbinfoVec)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // TinyExpr parser setup (as in setperiodicfaces)
    double x = 0.0, y = 0.0, z = 0.0;
    te_parser tep;
    tep.set_variables_and_functions({{"x", &x}, {"y", &y}, {"z", &z}});

    // Iterate over each periodic pair
    for (int ip = 0; ip < nprd; ++ip) {
        // Boundary-expression indices (0-based) for the two sides
        const int bc1 = prd_f1[ip] - 1;
        const int bc2 = prd_f2[ip] - 1;

        // ----------------------------------------------
        // 1. Collect local periodic faces for this pair
        // ----------------------------------------------
        // meta: [rank, elemLocal, elemGlobal, lf] per face
        std::vector<int>    meta1Local;
        std::vector<double> bary1Local;  // [n1loc * dim]

        std::vector<int>    meta2Local;
        std::vector<double> bary2Local;  // [n2loc * dim]

        meta1Local.reserve(nbfaces * 4);
        meta2Local.reserve(nbfaces * 4);
        bary1Local.reserve(nbfaces * dim);
        bary2Local.reserve(nbfaces * dim);

        for (int e = 0; e < ne; ++e) {
            int eGlob = elemGlobalID[e];

            for (int lf = 0; lf < nfe; ++lf) {
                int pos = e * nfe + lf;
                int val = t2t[pos];

                if (val >= 0) {
                    // already has a (non-boundary) neighbor
                    continue;
                }

                // t2t = -(k+1) encodes boundary label k
                int k = -val - 1;
                if (k != bc1 && k != bc2) {
                    continue; // not part of this periodic pair
                }

                // Compute face barycenter
                double xb[3] = {0.0, 0.0, 0.0};
                for (int i = 0; i < nvf; ++i) {
                    int ln_face = localfaces[lf * nvf + i];  // local node in element
                    int node    = t[e * nve + ln_face];      // local node index in mesh
                    const double* xp = &p[dim * node];
                    for (int d = 0; d < dim; ++d) {
                        xb[d] += xp[d];
                    }
                }
                double inv = 1.0 / static_cast<double>(nvf);
                for (int d = 0; d < dim; ++d) {
                    xb[d] *= inv;
                }

                if (k == bc1) {
                    meta1Local.insert(meta1Local.end(), {rank, e, eGlob, lf});
                    // meta1Local.push_back(rank);
                    // meta1Local.push_back(e);
                    // meta1Local.push_back(eGlob);
                    // meta1Local.push_back(lf);
                    for (int d = 0; d < dim; ++d)
                        bary1Local.push_back(xb[d]);
                } else { // k == bc2
                    meta2Local.insert(meta2Local.end(), {rank, e, eGlob, lf});
                    // meta2Local.push_back(rank);
                    // meta2Local.push_back(e);
                    // meta2Local.push_back(eGlob);
                    // meta2Local.push_back(lf);
                    for (int d = 0; d < dim; ++d)
                        bary2Local.push_back(xb[d]);
                }
            }
        }

        int n1loc = static_cast<int>(meta1Local.size() / 4);
        int n2loc = static_cast<int>(meta2Local.size() / 4);

        // ---------------------------------------------
        // 2. Evaluate mapping expressions at barycenters
        // ---------------------------------------------
        std::vector<double> q1Local(n1loc * ncomp, 0.0);
        std::vector<double> q2Local(n2loc * ncomp, 0.0);

        // Side 1: expr1
        for (int c = 0; c < ncomp; ++c) {
            auto result = tep.evaluate(expr1[ip * ncomp + c]);
            if (!tep.success()) {
                std::cout << "\t " << std::setfill(' ')
                          << std::setw(tep.get_last_error_position()) << '^'
                          << "\tError near here in periodic expr1\n";
                error("TinyExpr Failure in periodic expr1");
            }

            for (int j = 0; j < n1loc; ++j) {
                const double* xb = &bary1Local[j * dim];
                x = (dim > 0) ? xb[0] : 0.0;
                y = (dim > 1) ? xb[1] : 0.0;
                z = (dim > 2) ? xb[2] : 0.0;
                q1Local[j * ncomp + c] = tep.evaluate(expr1[ip * ncomp + c]);
            }
        }

        // Side 2: expr2
        for (int c = 0; c < ncomp; ++c) {
            auto result = tep.evaluate(expr2[ip * ncomp + c]);
            if (!tep.success()) {
                std::cout << "\t " << std::setfill(' ')
                          << std::setw(tep.get_last_error_position()) << '^'
                          << "\tError near here in periodic expr2\n";
                error("TinyExpr Failure in periodic expr2");
            }

            for (int j = 0; j < n2loc; ++j) {
                const double* xb = &bary2Local[j * dim];
                x = (dim > 0) ? xb[0] : 0.0;
                y = (dim > 1) ? xb[1] : 0.0;
                z = (dim > 2) ? xb[2] : 0.0;
                q2Local[j * ncomp + c] = tep.evaluate(expr2[ip * ncomp + c]);
            }
        }

        // ---------------------------------------
        // 3. MPI Allgather periodic faces for pair
        // ---------------------------------------
        std::vector<int> counts1(size), counts2(size);
        int n1loc_faces = n1loc;
        int n2loc_faces = n2loc;

        MPI_Allgather(&n1loc_faces, 1, MPI_INT, counts1.data(), 1, MPI_INT, comm);
        MPI_Allgather(&n2loc_faces, 1, MPI_INT, counts2.data(), 1, MPI_INT, comm);

        std::vector<int> disp1(size), disp2(size);
        int n1glob = 0, n2glob = 0;
        for (int r = 0; r < size; ++r) {
            disp1[r] = n1glob;
            disp2[r] = n2glob;
            n1glob += counts1[r];
            n2glob += counts2[r];
        }

        std::vector<int>    meta1_g(4 * n1glob), meta2_g(4 * n2glob);
        std::vector<double> q1_g(ncomp * n1glob), q2_g(ncomp * n2glob);

        // Allgatherv meta1
        {
            std::vector<int> metaCounts1(size), metaDisp1(size);
            for (int r = 0; r < size; ++r) {
                metaCounts1[r] = counts1[r] * 4;
                metaDisp1[r]   = disp1[r]   * 4;
            }
            MPI_Allgatherv(meta1Local.data(), n1loc * 4, MPI_INT,
                           meta1_g.data(), metaCounts1.data(), metaDisp1.data(), MPI_INT,
                           comm);
        }

        // Allgatherv meta2
        {
            std::vector<int> metaCounts2(size), metaDisp2(size);
            for (int r = 0; r < size; ++r) {
                metaCounts2[r] = counts2[r] * 4;
                metaDisp2[r]   = disp2[r]   * 4;
            }
            MPI_Allgatherv(meta2Local.data(), n2loc * 4, MPI_INT,
                           meta2_g.data(), metaCounts2.data(), metaDisp2.data(), MPI_INT,
                           comm);
        }

        // Allgatherv q1
        {
            std::vector<int> qCounts1(size), qDisp1(size);
            for (int r = 0; r < size; ++r) {
                qCounts1[r] = counts1[r] * ncomp;
                qDisp1[r]   = disp1[r]   * ncomp;
            }
            MPI_Allgatherv(q1Local.data(), n1loc * ncomp, MPI_DOUBLE,
                           q1_g.data(), qCounts1.data(), qDisp1.data(), MPI_DOUBLE,
                           comm);
        }

        // Allgatherv q2
        {
            std::vector<int> qCounts2(size), qDisp2(size);
            for (int r = 0; r < size; ++r) {
                qCounts2[r] = counts2[r] * ncomp;
                qDisp2[r]   = disp2[r]   * ncomp;
            }
            MPI_Allgatherv(q2Local.data(), n2loc * ncomp, MPI_DOUBLE,
                           q2_g.data(), qCounts2.data(), qDisp2.data(), MPI_DOUBLE,
                           comm);
        }

        if (n1glob == 0 || n2glob == 0) {
            // no periodic faces of this pair on this mesh
            continue;
        }

        // ---------------------------------------
        // 4. Match faces using xiny on q1 vs q2
        // ---------------------------------------
        std::vector<int> in(n1glob, -1);
        xiny<double>(in.data(), q1_g.data(), q2_g.data(),
                     n1glob, n2glob, ncomp);

        // ---------------------------------------
        // 5. Update t2t and nbinfoVec
        // ---------------------------------------
        for (int j = 0; j < n1glob; ++j) {
            int idx2 = in[j];
            if (idx2 < 0) continue; // no match

            // side 1
            int r1     = meta1_g[4*j + 0];
            int eLoc1  = meta1_g[4*j + 1];
            int eGlob1 = meta1_g[4*j + 2];
            int lf1    = meta1_g[4*j + 3];

            // side 2
            int r2     = meta2_g[4*idx2 + 0];
            int eLoc2  = meta2_g[4*idx2 + 1];
            int eGlob2 = meta2_g[4*idx2 + 2];
            int lf2    = meta2_g[4*idx2 + 3];

            // Update side 1 on its owning rank:
            if (r1 == rank) {
                int pos1 = eLoc1 * nfe + lf1;
                t2t[pos1] = eGlob2;  // periodic neighbor is global element ID

                if (r1 != r2) {
                    nbinfoVec.insert(nbinfoVec.end(), {eLoc1, lf1, rank, eGlob2, lf2, r2});
                    // // cross-rank neighbor: record in nbinfoVec
                    // nbinfoVec.push_back(eLoc1);   // self element local
                    // nbinfoVec.push_back(lf1);     // self local face
                    // nbinfoVec.push_back(rank);    // self rank
                    // nbinfoVec.push_back(eGlob2);  // neighbor element global
                    // nbinfoVec.push_back(lf2);     // neighbor local face (on its rank)
                    // nbinfoVec.push_back(r2);      // neighbor rank
                }
            }

            // Update side 2 on its owning rank:
            if (r2 == rank) {
                int pos2 = eLoc2 * nfe + lf2;
                t2t[pos2] = eGlob1;

                if (r2 != r1) {
                    nbinfoVec.insert(nbinfoVec.end(), {eLoc2, lf2, rank, eGlob1, lf1, r1});
                    // nbinfoVec.push_back(eLoc2);
                    // nbinfoVec.push_back(lf2);
                    // nbinfoVec.push_back(rank);
                    // nbinfoVec.push_back(eGlob1);
                    // nbinfoVec.push_back(lf1);
                    // nbinfoVec.push_back(r1);
                }
            }
        } // j
    } // ip
}

void computeElementToGlobalNodeMap(Mesh& mesh, const vector<int> &elempart)
{
    const int nve = mesh.nve;
    const int ne  = mesh.ne;

    // Allocate output
    mesh.tg.resize(static_cast<size_t>(nve) * ne);

    // Loop over elements
    for (int e = 0; e < ne; ++e) {

        const int offset = elempart[e] * nve;

        // Map each local node to its global node ID
        for (int i = 0; i < nve; ++i) {
            int localNode = mesh.t[offset + i];
            mesh.tg[i + nve*e] = mesh.nodeGlobalID[localNode];
        }
    }
}

void sortWithReordering(std::vector<int>& a, std::vector<int>& b)
{
    // Ensure same size
    if (a.size() != b.size()) {
        throw std::runtime_error("a and b must have the same size");
    }

    // Create index array
    std::vector<std::size_t> idx(a.size());
    for (std::size_t i = 0; i < idx.size(); ++i) {
        idx[i] = i;
    }

    // Sort indices based on values in a
    std::sort(idx.begin(), idx.end(),
              [&](std::size_t i, std::size_t j) { return a[i] < a[j]; });

    // Create sorted arrays
    std::vector<int> a_sorted(a.size());
    std::vector<int> b_sorted(b.size());

    for (std::size_t k = 0; k < idx.size(); ++k) {
        a_sorted[k] = a[idx[k]];
        b_sorted[k] = b[idx[k]];
    }

    // Replace original vectors
    a = std::move(a_sorted);
    b = std::move(b_sorted);
}

// nbinfo layout per record: 6 ints
void classifyElementsWithE2EAndNbinfo(const int* e2e,       // [mesh.ne * nfe], global IDs or -1
                                      const int* elemGlobalID,
                                      int        nfe,
                                      int        ne_local,
                                      const int* nbinfo,    // may be nullptr if no remote faces
                                      int        nbinfoSize,
                                      ElementClassification& out, int rank)
{
    out.interiorLocal.clear();
    out.boundaryLocal.clear();
    out.interfaceLocal.clear();
    out.interiorGlobal.clear();
    out.boundaryGlobal.clear();
    out.interfaceGlobal.clear();
    out.neighborElemGlobal.clear();
    out.neighborElemRank.clear();

    if (ne_local == 0)
        return;

    // --------------------------------------------------------------
    // 1. Flags: does element have a physical boundary face?
    //           does element have a remote (off-rank) face?
    // --------------------------------------------------------------
    std::vector<char> hasBoundaryFace(ne_local, 0);
    std::vector<char> hasRemoteFace(ne_local,   0);

    // Boundary faces: e2e[e*nfe + lf] <= -1 means physical boundary
    for (int e = 0; e < ne_local; ++e) {
        for (int lf = 0; lf < nfe; ++lf) {
            int ngh = e2e[e * nfe + lf];
            if (ngh <= -1) {
                hasBoundaryFace[e] = 1;
            }
        }
    }

    // Remote faces & neighbor elements come from nbinfo
    // nbinfo record size:
    const int nbRecordSize = 6;

    // To deduplicate remote neighbor elements
    struct NeighborKey {
        int globalID;
        int rank;
        bool operator==(const NeighborKey& o) const noexcept {
            return globalID == o.globalID && rank == o.rank;
        }
    };
    struct NeighborKeyHash {
        std::size_t operator()(const NeighborKey& k) const noexcept {
            std::size_t h1 = std::hash<int>{}(k.globalID);
            std::size_t h2 = std::hash<int>{}(k.rank);
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
    std::unordered_set<NeighborKey, NeighborKeyHash> neighborSet;

    for (int r = 0; r < nbinfoSize; ++r) {
        const int base = r * nbRecordSize;

        int eLoc      = nbinfo[base + 0];
        int lf        = nbinfo[base + 1];
        int selfRank  = nbinfo[base + 2]; (void)selfRank; // not used here
        int eglonb    = nbinfo[base + 3];
        int lfnb      = nbinfo[base + 4]; (void)lfnb;     // neighbor face index, not needed for classification
        int nghRank   = nbinfo[base + 5];

        // Mark this element as having a remote face
        if (eLoc < 0 || eLoc >= ne_local) {
            // optionally error() here if you enforce strict checks
            continue;
        }
        hasRemoteFace[eLoc] = 1;

        // Record this remote neighbor element (deduplicated)
        NeighborKey key{eglonb, nghRank};
        if (neighborSet.insert(key).second) {
            out.neighborElemGlobal.push_back(eglonb);
            out.neighborElemRank.push_back(nghRank);
        }

        // We ignore the neighbor's e2e row in nbinfo[base + 6 .. base + 6+nfe-1]
        // for pure classification purposes, but it's available if needed.
    }

    sortWithReordering(out.neighborElemGlobal, out.neighborElemRank);

    // --------------------------------------------------------------
    // 2. Classify each local element
    //
    // Convention (disjoint sets):
    //   - interface elements: hasRemoteFace[e] == 1
    //   - boundary elements: !hasRemoteFace[e] && hasBoundaryFace[e] == 1
    //   - interior elements: !hasRemoteFace[e] && !hasBoundaryFace[e]
    //
    // If you want elements that touch both remote and boundary faces to be
    // in *both* sets, you can change the logic accordingly.
    // --------------------------------------------------------------
    for (int e = 0; e < ne_local; ++e) {
        int g = elemGlobalID[e];

        if (hasRemoteFace[e]) {
            out.interfaceLocal.push_back(e);
            out.interfaceGlobal.push_back(g);
        } else if (hasBoundaryFace[e]) {
            out.boundaryLocal.push_back(e);
            out.boundaryGlobal.push_back(g);
        } else {
            out.interiorLocal.push_back(e);
            out.interiorGlobal.push_back(g);
        }
    }

    if (rank == 0) {
        std::cout << "Finished classifying elements on each subdomain" << std::endl;
    }      
}


void buildElempartFromClassification(ElementClassification& cls, DMD& dmd, int rank)
{
    // -------------------------------------------------------------
    // 1. Clear previous data
    // -------------------------------------------------------------
    dmd.elempart.clear();
    dmd.elempartpts.clear();
    dmd.elempart_local.clear();

    // Counts
    const int nInterior = cls.interiorGlobal.size();
    const int nBoundary = cls.boundaryGlobal.size();
    const int nInterface = cls.interfaceGlobal.size();
    const int nNeighbor = cls.neighborElemGlobal.size();
    const int nOuter = cls.outerElemGlobal.size();

    // -------------------------------------------------------------
    // 2. Concatenate into elempart in the desired order
    //
    // interiorGlobal   – global interior elements
    // boundaryGlobal   – global boundary elements
    // interfaceGlobal  – global interface (shared) elements
    // neighborElemGlobal – off-rank neighbor elements
    // -------------------------------------------------------------
    dmd.elempart.reserve(nInterior + nBoundary + nInterface + nNeighbor + nOuter);

    //  boundary
    dmd.elempart.insert(dmd.elempart.end(),
                        cls.boundaryGlobal.begin(),
                        cls.boundaryGlobal.end());

    //  interior
    dmd.elempart.insert(dmd.elempart.end(),
                        cls.interiorGlobal.begin(),
                        cls.interiorGlobal.end());

    //  interface
    dmd.elempart.insert(dmd.elempart.end(),
                        cls.interfaceGlobal.begin(),
                        cls.interfaceGlobal.end());

    //  off-rank neighbor
    dmd.elempart.insert(dmd.elempart.end(),
                        cls.neighborElemGlobal.begin(),
                        cls.neighborElemGlobal.end());

    cls.neighborElemLocal.resize(nNeighbor);
    const int neighborOffset = nBoundary + nInterior + nInterface;
    for (int i = 0; i < nNeighbor; ++i) {
        cls.neighborElemLocal[i] = neighborOffset + i;
    }

    // -------------------------------------------------------------
    // 2b. Build elempart_local in the same order
    //
    // elempart_local holds local element indices corresponding to
    // the global IDs in elempart.
    // -------------------------------------------------------------
    dmd.elempart_local.reserve(nInterior + nBoundary + nInterface + nNeighbor + nOuter);

    // boundary (local)
    dmd.elempart_local.insert(dmd.elempart_local.end(),
                              cls.boundaryLocal.begin(),
                              cls.boundaryLocal.end());

    // interior (local)
    dmd.elempart_local.insert(dmd.elempart_local.end(),
                              cls.interiorLocal.begin(),
                              cls.interiorLocal.end());

    // interface (local)
    dmd.elempart_local.insert(dmd.elempart_local.end(),
                              cls.interfaceLocal.begin(),
                              cls.interfaceLocal.end());

    // off-rank neighbor (local indices of ghost / halo elements)
    dmd.elempart_local.insert(dmd.elempart_local.end(),
                              cls.neighborElemLocal.begin(),
                              cls.neighborElemLocal.end());  

    if (nOuter > 0) {
      dmd.elempart.insert(dmd.elempart.end(),
                        cls.outerElemGlobal.begin(),
                        cls.outerElemGlobal.end());      

      cls.outerElemLocal.resize(nOuter);
      const int outerOffset = nBoundary + nInterior + nInterface + nNeighbor;
      for (int i = 0; i < nOuter; ++i) cls.outerElemLocal[i] = outerOffset + i;
            
      dmd.elempart_local.insert(dmd.elempart_local.end(),
                          cls.outerElemLocal.begin(),
                          cls.outerElemLocal.end());        
    }
  
    // -------------------------------------------------------------
    // 3. Fill elempartpts
    //
    // You defined elempartpts as:
    //     [interior, interface, exterior]
    //
    // but boundary elements belong to neither interior nor interface
    // -------------------------------------------------------------
    dmd.elempartpts.resize(3);
    dmd.elempartpts[0] = nInterior + nBoundary;// interior + boundary
    dmd.elempartpts[1] = nInterface;           // interface
    dmd.elempartpts[2] = nNeighbor;            // exterior

    if (nOuter > 0) dmd.elempartpts.push_back(nOuter);

    // Optional: intepartpts if you want a 4-way split (interior, boundary, interface, exterior)
    if (!dmd.intepartpts.empty()) {
        dmd.intepartpts.resize(4);
        dmd.intepartpts[0] = nBoundary;
        dmd.intepartpts[1] = nInterior;
        dmd.intepartpts[2] = nInterface;
        dmd.intepartpts[3] = nNeighbor;
        if (nOuter > 0) dmd.intepartpts.push_back(nOuter);
    }

    if (rank == 0) {
        std::cout << "Finished partitioning elements on each subdomain" << std::endl;
    }        
}

void buildElem2CpuFromClassification(const ElementClassification& cls,
                                     DMD&                        dmd,
                                     MPI_Comm                    comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    (void)size;

    const int nInterior  = static_cast<int>(cls.interiorGlobal.size());
    const int nBoundary  = static_cast<int>(cls.boundaryGlobal.size());
    const int nInterface = static_cast<int>(cls.interfaceGlobal.size());
    const int nNeighbor  = static_cast<int>(cls.neighborElemGlobal.size());
    const int nOuter  = static_cast<int>(cls.outerElemGlobal.size());

    const int expectedSize = nInterior + nBoundary + nInterface + nNeighbor + nOuter;
    const int elempartSize = static_cast<int>(dmd.elempart.size());

    // Sanity check: layout of dmd.elempart must match classification
    if (elempartSize != expectedSize) {
        std::cerr << "buildElem2CpuFromClassification: size mismatch:\n"
                  << "  dmd.elempart.size() = " << elempartSize << "\n"
                  << "  expected            = " << expectedSize << "\n";
        // You can replace with your own error() routine if you prefer.
        // error("buildElem2CpuFromClassification: elempart mismatch");
    }

    dmd.elem2cpu.resize(elempartSize);

    int pos = 0;

    // 1) boundaryGlobal → also owned by this rank
    for (int i = 0; i < nBoundary; ++i) {
        dmd.elem2cpu[pos++] = rank;
    }
  
    // 2) interiorGlobal → owned by this rank
    for (int i = 0; i < nInterior; ++i) {
        dmd.elem2cpu[pos++] = rank;
    }

    // 3) interfaceGlobal → also owned by this rank
    for (int i = 0; i < nInterface; ++i) {
        dmd.elem2cpu[pos++] = rank;
    }

    // 4) neighborElemGlobal → owned by their respective neighbor ranks
    //    Use cls.neighborElemRank (same ordering as neighborElemGlobal).
    for (int i = 0; i < nNeighbor; ++i) {
        int ownerRank = cls.neighborElemRank[i];
        // optional sanity check
        if (ownerRank < 0 || ownerRank >= size) {
            std::cerr << "buildElem2CpuFromClassification: invalid owner rank "
                      << ownerRank << " for neighbor element " << i << "\n";
        }
        dmd.elem2cpu[pos++] = ownerRank;
    }

    // 5) outerElemGlobal → owned by their respective outer ranks
    //    Use cls.outerElemRank (same ordering as outerElemGlobal).
    for (int i = 0; i < nOuter; ++i) {
        int ownerRank = cls.outerElemRank[i];
        // optional sanity check
        if (ownerRank < 0 || ownerRank >= size) {
            std::cerr << "buildElem2CpuFromClassification: invalid owner rank "
                      << ownerRank << " for neighbor element " << i << "\n";
        }
        dmd.elem2cpu[pos++] = ownerRank;
    }
  
    // Final consistency check
    assert(pos == elempartSize);

    if (rank == 0) {
        std::cout << "Finished building Elem2Rank on each subdomain" << std::endl;
    }          
}

void buildElemRecv(DMD& dmd, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    (void)size;

    dmd.elemrecv.clear();
    dmd.elemrecvpts.clear();
    dmd.nbsd.clear();

    if (dmd.elempart.empty() || dmd.elem2cpu.empty() || dmd.elempartpts.size() < 3) {
        return; // nothing to do
    }

    const int nInterior  = dmd.elempartpts[0];
    const int nInterface = dmd.elempartpts[1];
    int nExterior  = dmd.elempartpts[2];
    if (dmd.elempartpts.size() > 3) nExterior += dmd.elempartpts[3];

    const int nPart     = static_cast<int>(dmd.elempart.size());
    const int expected  = nInterior + nInterface + nExterior;

    if (nPart != expected || static_cast<int>(dmd.elem2cpu.size()) != nPart) {
        std::cerr << "buildElemRecv: inconsistent elempart / elem2cpu / elempartpts\n";
        return;
    }

    // Reserve memory based on maximum possible sizes
    dmd.elemrecv.reserve(nExterior);     // only exterior elements produce ghosts
    dmd.nbsd.reserve(nExterior);         // worst-case: each ghost from a different rank
    dmd.elemrecvpts.reserve(nExterior);  // one entry per neighbor
  
    // Range of ghost (exterior) elements in elempart
    const int extStart = nInterior + nInterface;
    const int extEnd   = extStart + nExterior; // == nPart

    // --------------------------------------------------------
    // 1. Build elemrecv rows: [owner, recv_local_idx, globalID]
    // --------------------------------------------------------
    for (int k = extStart; k < extEnd; ++k) {
        int owner    = dmd.elem2cpu[k];     // owning rank of this element
        int globalID = dmd.elempart[k];     // global element ID

        if (owner == rank) {
            // Shouldn't normally happen for "exterior" elements,
            // but ignore quietly if it does.
            continue;
        }
        if (owner < 0) {
            std::cerr << "buildElemRecv: negative owner rank at elempart index "
                      << k << "\n";
            continue;
        }

        dmd.elemrecv.push_back({owner, k, globalID});
    }

    if (dmd.elemrecv.empty()) {
        // no ghosts → no neighbors
        return;
    }
    

    //if (rank==0) printElemRecv(dmd.elemrecv);
  
    // --------------------------------------------------------
    // 2. (Optional but usually helpful) group rows by owner rank
    // --------------------------------------------------------
    std::sort(dmd.elemrecv.begin(), dmd.elemrecv.end());
    // std::sort(dmd.elemrecv.begin(), dmd.elemrecv.end(),
    //           [](const std::vector<int>& a, const std::vector<int>& b) {
    //               // row = [owner, recv_local_idx, globalID]
    //               return a[0] < b[0]; // sort by owner
    //           });
    // std::sort(dmd.elemrecv.begin(), dmd.elemrecv.end(),
    //     [](const std::array<int, 3>& a, const std::array<int, 3>& b) {
    //         return a[0] < b[0];   // sort by owner
    //     });
  
    // --------------------------------------------------------
    // 3. Build nbsd (list of neighbor ranks) and elemrecvpts
    //    (number of elements received from each neighbor)
    // --------------------------------------------------------
    int currentOwner = -1;
    int countForOwner = 0;

    for (std::size_t i = 0; i < dmd.elemrecv.size(); ++i) {
        int owner = dmd.elemrecv[i][0];

        if (owner != currentOwner) {
            // flush previous owner count
            if (currentOwner != -1) {
                dmd.nbsd.push_back(currentOwner);
                dmd.elemrecvpts.push_back(countForOwner);
            }
            currentOwner = owner;
            countForOwner = 1;
        } else {
            ++countForOwner;
        }
    }

    // flush last owner
    if (currentOwner != -1) {
        dmd.nbsd.push_back(currentOwner);
        dmd.elemrecvpts.push_back(countForOwner);
    }

    if (rank == 0) {
        std::cout << "Finished building ElemRecv on each subdomain" << std::endl;
    }        
  
    // At this point:
    //   - dmd.elemrecv is sorted by owner; each row = [owner, recv_local_idx, globalID]
    //   - dmd.nbsd[i] is neighbor rank i
    //   - dmd.elemrecvpts[i] is #ghost elements received from dmd.nbsd[i]
}

void buildElemsend(const Mesh& mesh, DMD& dmd, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    dmd.elemsend.clear();
    dmd.elemsendpts.clear();

    const int ne_local   = mesh.ne;
    const int nRecvRows  = static_cast<int>(dmd.elemrecv.size());
    const auto& nbsd     = dmd.nbsd;  // neighbor list
    const int nNeighbors = static_cast<int>(nbsd.size());

    if (ne_local == 0 && nRecvRows == 0)
        return;

    // Fast lookup: is this rank a neighbor?
    std::vector<char> isNeighbor(size, 0);
    for (int nbr : nbsd) {
        if (nbr >= 0 && nbr < size)
            isNeighbor[nbr] = 1;
    }

    // ------------------------------------------------------------
    // 1. Build globalID → localID map for *owned* elements
    // ------------------------------------------------------------
    std::unordered_map<int,int> gid2loc;
    gid2loc.reserve(static_cast<size_t>(ne_local) * 2);

    for (int eLoc = 0; eLoc < ne_local; ++eLoc) {
      gid2loc[ dmd.elempart[eLoc] ] = eLoc;
    }

    // ------------------------------------------------------------
    // 2. Count how many requests we send to each rank
    //    elemrecv row: [owner_rank, recv_local_idx, globalID]
    // ------------------------------------------------------------
    std::vector<int> sendCounts(size, 0);

    for (const auto& row : dmd.elemrecv)
    {
        const int owner = row[0];

        if (owner < 0 || owner >= size) continue;
        if (!isNeighbor[owner])         continue;
        if (owner == rank)              continue; // should not happen, but safe

        sendCounts[owner] += 1;
    }

    // ------------------------------------------------------------
    // 3. Nonblocking exchange of request counts
    // ------------------------------------------------------------
    constexpr int TAG_COUNTS = 9001;
    std::vector<int> recvCounts(size, 0);

    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * nNeighbors);

    // Irecv from neighbors
    for (int nbr : nbsd) {
        MPI_Request rq;
        MPI_Irecv(&recvCounts[nbr], 1, MPI_INT, nbr, TAG_COUNTS, comm, &rq);
        reqs.push_back(rq);
    }

    // Isend to neighbors
    for (int nbr : nbsd) {
        MPI_Request rq;
        MPI_Isend(&sendCounts[nbr], 1, MPI_INT, nbr, TAG_COUNTS, comm, &rq);
        reqs.push_back(rq);
    }

    if (!reqs.empty())
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);

    // ------------------------------------------------------------
    // 4. Allocate send/recv buffers
    // ------------------------------------------------------------
    std::vector<std::vector<int>> sendReqBuf(size);
    std::vector<std::vector<int>> recvReqBuf(size);

    for (int nbr : nbsd) {
        sendReqBuf[nbr].resize(sendCounts[nbr]);
        recvReqBuf[nbr].resize(recvCounts[nbr]);
    }

    // ------------------------------------------------------------
    // 5. Pack the request buffers (global element IDs)
    // ------------------------------------------------------------
    std::vector<int> offset(size, 0);

    for (const auto& row : dmd.elemrecv)
    {
        const int owner = row[0];
        const int gID   = row[2];

        if (owner < 0 || owner >= size) continue;
        if (!isNeighbor[owner])         continue;

        const int pos = offset[owner]++;
        sendReqBuf[owner][pos] = gID;
    }

    // ------------------------------------------------------------
    // 6. Nonblocking exchange of request data
    // ------------------------------------------------------------
    reqs.clear();

    constexpr int TAG_REQS = 9002;

    // Irecv
    for (int nbr : nbsd) {
        const int count = recvCounts[nbr];
        if (count == 0) continue;

        MPI_Request rq;
        MPI_Irecv(recvReqBuf[nbr].data(),
                  count, MPI_INT,
                  nbr, TAG_REQS, comm, &rq);
        reqs.push_back(rq);
    }

    // Isend
    for (int nbr : nbsd) {
        const int count = sendCounts[nbr];
        if (count == 0) continue;

        MPI_Request rq;
        MPI_Isend(sendReqBuf[nbr].data(),
                  count, MPI_INT,
                  nbr, TAG_REQS, comm, &rq);
        reqs.push_back(rq);
    }

    if (!reqs.empty())
        MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);

    // ------------------------------------------------------------
    // 7. Reserve elemsend and build it from received requests
    //    Also compute per-neighbor counts, then build cumulative elemsendpts.
    // ------------------------------------------------------------
    std::size_t totalSend = 0;
    for (int nbr : nbsd) {
        totalSend += static_cast<std::size_t>(recvCounts[nbr]);
    }
    dmd.elemsend.reserve(totalSend);

    std::vector<int> sendPerNeighbor(nNeighbors, 0); // count per nbsd[i]

    for (int iNbr = 0; iNbr < nNeighbors; ++iNbr)
    {
        const int nbr  = nbsd[iNbr];
        const int nReq = recvCounts[nbr];
        if (nReq == 0) continue;

        const auto& reqsFromNbr = recvReqBuf[nbr];

        for (int k = 0; k < nReq; ++k)
        {
            const int gID = reqsFromNbr[k];

            auto it = gid2loc.find(gID);
            if (it == gid2loc.end()) {
                std::cerr << "ERROR: owner rank " << rank
                          << " does not own requested elem " << gID
                          << " from rank " << nbr << "\n";
                continue;
            }

            const int localIdx = it->second;

            // elemsend row: [receiver_rank, local_idx, globalID]
            dmd.elemsend.push_back({nbr, localIdx, gID});
            ++sendPerNeighbor[iNbr];
        }
    }

    // ------------------------------------------------------------
    // 8. Build cumulative elemsendpts (size = nNeighbors + 1)
    //    elemsendpts[i] = total #elems sent to neighbors 0..i-1
    // ------------------------------------------------------------
    dmd.elemsendpts.resize(nNeighbors);
    for (int i = 0; i < nNeighbors; ++i) {
        dmd.elemsendpts[i] = sendPerNeighbor[i];
    }
    // dmd.elemsendpts.resize(nNeighbors + 1);
    // dmd.elemsendpts[0] = 0;
    // for (int i = 0; i < nNeighbors; ++i) {
    //     dmd.elemsendpts[i + 1] = dmd.elemsendpts[i] + sendPerNeighbor[i];
    // }

    if (rank == 0) {
        std::cout << "Finished building ElemSend on each subdomain" << std::endl;
    }
}

// Helper: map C++ type to MPI_Datatype
template<typename T>
MPI_Datatype mpi_type();

template<>
inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }

template<>
inline MPI_Datatype mpi_type<float>()  { return MPI_FLOAT; }

template<>
inline MPI_Datatype mpi_type<int>()    { return MPI_INT; }

template<typename T>
void sendrecvdata(
    MPI_Comm comm,
    const std::vector<int>   &nbsd,       // neighbor ranks (size = nnbsd)
    const std::vector<int> &elemsendpts,  // #elements to send to each neighbor
    const std::vector<int> &elemrecvpts,  // #elements to recv from each neighbor
    const std::vector<int> &elemsend,     // list of local element IDs to send
    const std::vector<int> &elemrecv,     // list of local element IDs to receive
    std::vector<T>         &senddata,     // full source data (all elements)
    std::vector<T>         &recvdata,     // full destination data (all elements)
    int bsz                               // #entries per element (block size)
)
{
    const int nelemsend = static_cast<int>(elemsend.size());
    const int nelemrecv = static_cast<int>(elemrecv.size());
    const int   nnbsd     = static_cast<int>(nbsd.size());

    // -------------------------------------------------------------------------
    // 0. Local communication buffers
    // -------------------------------------------------------------------------
    std::vector<T> buffsend(static_cast<std::size_t>(nelemsend) * bsz);
    std::vector<T> buffrecv(static_cast<std::size_t>(nelemrecv) * bsz);

    // -------------------------------------------------------------------------
    // 1. Pack send buffer from senddata according to elemsend
    // -------------------------------------------------------------------------
    select_columns(
        buffsend.data(),        // packed send buffer (contiguous)
        senddata.data(),        // global/source data
        elemsend.data(),        // element indices to extract
        bsz,
        nelemsend
    );

    // -------------------------------------------------------------------------
    // 2. Count messages
    // -------------------------------------------------------------------------
    int nsends = 0, nrecvs = 0;
    for (int n = 0; n < nnbsd; ++n) {
        if (elemsendpts[n] * static_cast<int>(bsz) > 0) ++nsends;
        if (elemrecvpts[n] * static_cast<int>(bsz) > 0) ++nrecvs;
    }

    const int totalReq = nsends + nrecvs;
    std::vector<MPI_Request> requests(totalReq);
    std::vector<MPI_Status>  statuses(totalReq);

    int request_counter = 0;
    const MPI_Datatype mpiT = mpi_type<T>();

    // -------------------------------------------------------------------------
    // 3. Non-blocking sends
    // -------------------------------------------------------------------------
    int psend = 0;
    for (int n = 0; n < nnbsd; ++n) {
        const int   neighbor = nbsd[n];
        const int nsend    = elemsendpts[n] * static_cast<int>(bsz);

        if (nsend > 0) {
            MPI_Isend(
                &buffsend[static_cast<std::size_t>(psend)],
                static_cast<int>(nsend),
                mpiT,
                neighbor,
                0,                    // tag
                comm,
                &requests[request_counter++]
            );
            psend += nsend;
        }
    }

    // -------------------------------------------------------------------------
    // 4. Non-blocking receives
    // -------------------------------------------------------------------------
    int precv = 0;
    for (int n = 0; n < nnbsd; ++n) {
        const int   neighbor = nbsd[n];
        const int nrecv    = elemrecvpts[n] * static_cast<int>(bsz);

        if (nrecv > 0) {
            MPI_Irecv(
                &buffrecv[static_cast<std::size_t>(precv)],
                static_cast<int>(nrecv),
                mpiT,
                neighbor,
                0,                    // tag
                comm,
                &requests[request_counter++]
            );
            precv += nrecv;
        }
    }

    // -------------------------------------------------------------------------
    // 5. Wait for all communication to complete
    // -------------------------------------------------------------------------
    MPI_Waitall(request_counter, requests.data(), statuses.data());

    // -------------------------------------------------------------------------
    // 6. Unpack recv buffer back into recvdata according to elemrecv
    // -------------------------------------------------------------------------
    insert_columns(
        recvdata.data(),        // destination/global array
        buffrecv.data(),        // contiguous received data
        elemrecv.data(),        // element indices we received
        bsz,
        nelemrecv
    );
}

// Step 1: collect global neighbor IDs for interface elements
//
// ifaceNeighGlobal has size interfaceLocal.size() * nfe, stored as
//   ifaceNeighGlobal[k * nfe + lf] = global neighbor ID for
//   element eLoc = cls.interfaceLocal[k], face lf (or -1 if boundary).
//
void collectInterfaceNeighborGlobals(const int*                e2e,
                                     int                       nfe,
                                     const ElementClassification& cls,
                                     std::vector<int>&         ifaceNeighGlobal)
{
    const std::vector<int>& interfaceLocal = cls.interfaceLocal;
    const int nInterface = static_cast<int>(interfaceLocal.size());

    ifaceNeighGlobal.resize(static_cast<std::size_t>(nInterface) * nfe);

    for (int k = 0; k < nInterface; ++k) {
        int eLoc = interfaceLocal[k];  // local element index on this rank
        for (int lf = 0; lf < nfe; ++lf) {
            ifaceNeighGlobal[k * nfe + lf] =
                e2e[eLoc * nfe + lf];   // global neighbor ID or -1
        }
    }
}

// Step 2: determine owning rank for each global element ID in ifaceNeighGlobal.
// Inputs:
//   mesh.elemGlobalID : local elements' global IDs (for building the owner map)
//   ifaceNeighGlobal  : flattened nInterface * nfe array from Step 1
//
// Output:
//   ifaceNeighOwner   : same size, ifaceNeighOwner[i] = rank that owns
//                       ifaceNeighGlobal[i], or -1 if invalid.
// ifaceNeighGlobal: size = interfaceLocal.size() * nfe (from step 1)
// ifaceNeighOwner:  same size, owner rank for each global ID (or -1 for invalid)
//
// Assumption: every non-negative entry in ifaceNeighGlobal belongs to either
//   mesh.elemGlobalID (local elements) or cls.neighborElemGlobal (1-hop neighbors).
//
void determineOwnersForInterfaceNeighborsLocal(const Mesh&                  mesh,
                                               const ElementClassification& cls,
                                               const std::vector<int>&      ifaceNeighGlobal,
                                               std::vector<int>&            ifaceNeighOwner,
                                               MPI_Comm                     comm)
{
    int rank;
    MPI_Comm_rank(comm, &rank);

    ifaceNeighOwner.resize(ifaceNeighGlobal.size());

    // Build a simple map globalID -> ownerRank using only local + neighbor info.
    std::unordered_map<int,int> ownerMap;
    ownerMap.reserve(mesh.ne + cls.neighborElemGlobal.size() * 2);

    // 1) Local elements: owned by this rank
    for (int e = 0; e < mesh.ne; ++e) {
        int gid = mesh.elemGlobalID[e];
        ownerMap[gid] = rank;
    }

    // 2) Neighbor elements: owned by neighbor ranks
    for (std::size_t i = 0; i < cls.neighborElemGlobal.size(); ++i) {
        int gid    = cls.neighborElemGlobal[i];
        int nghRnk = cls.neighborElemRank[i];
        ownerMap[gid] = nghRnk;
    }

    // 3) Fill ifaceNeighOwner
    for (std::size_t i = 0; i < ifaceNeighGlobal.size(); ++i) {
        int gid = ifaceNeighGlobal[i];

        if (gid < 0) {
            // -1 → physical boundary / no neighbor
            ifaceNeighOwner[i] = -1;
            continue;
        }

        auto it = ownerMap.find(gid);
        if (it == ownerMap.end()) {
            // In principle this should not happen under your assumption
            // You can choose to error() here if you want it strict.
            ifaceNeighOwner[i] = -1;
        } else {
            ifaceNeighOwner[i] = it->second;
        }
    }
}

// Exchange the relevant portions of ifaceNeighGlobal / ifaceNeighOwner
// between neighboring ranks, as indicated by nbinfo.
//
// Inputs:
//   mesh            : local mesh (only used for mesh.ne)
//   cls.interfaceLocal : list of local interface elements
//   ifaceNeighGlobal : size = interfaceLocal.size() * nfe
//   ifaceNeighOwner  : same size as ifaceNeighGlobal
//   nbinfo           : length = 6 * nbinfoSize, records:
//                      (eLoc, lf_self, self_rank, eglonb, lfnb, neighbor_rank)
//   nbinfoSize       : number of nbinfo records
//   nfe              : number of faces per element
//
// Outputs:
//   recvBuf          : flat int buffer of all records received from neighbors
//   recvCounts       : per-rank counts of ints in recvBuf
//   recvDispls       : per-rank displacements into recvBuf
//
// Record layout in recvBuf (per received face record):
//   [0]           = eglonb    (global element ID on *this* rank)
//   [1]           = lfnb      (local face index on *this* rank)
//   [2 .. 2+nfe-1]           = ifaceNeighGlobal_row_from_remote[nfe ints]
//   [2+nfe .. 2+2*nfe-1]     = ifaceNeighOwner_row_from_remote [nfe ints]
//
void exchangeInterfaceNeighborRows(const Mesh&                  mesh,
                                   const ElementClassification& cls,
                                   const std::vector<int>&      ifaceNeighGlobal,
                                   const std::vector<int>&      ifaceNeighOwner,
                                   const int*                   nbinfo,
                                   int                          nbinfoSize,
                                   int                          nfe,
                                   MPI_Comm                     comm,
                                   std::vector<int>&            recvBuf,
                                   std::vector<int>&            recvCounts,
                                   std::vector<int>&            recvDispls)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int ne_local = mesh.ne;
    const int nbRecordSize = 6;          // nbinfo record size
    const int ifaceRecordSize = 2 + 2*nfe; // eglonb, lfnb, G[nfe], Owner[nfe]

    if (ne_local == 0 || cls.interfaceLocal.empty() || nbinfoSize == 0) {
        recvBuf.clear();
        recvCounts.assign(size, 0);
        recvDispls.assign(size, 0);
        return;
    }

    // ------------------------------------------------------------------
    // 1. Build mapping: local element -> interface index
    //    so we can find the row in ifaceNeighGlobal/Owner for eLoc.
    // ------------------------------------------------------------------
    std::vector<int> local2iface(ne_local, -1);
    for (int k = 0; k < (int)cls.interfaceLocal.size(); ++k) {
        int eLoc = cls.interfaceLocal[k];
        if (eLoc >= 0 && eLoc < ne_local) {
            local2iface[eLoc] = k;
        }
    }

    // ------------------------------------------------------------------
    // 2. Count how many records we will send to each neighbor rank
    // ------------------------------------------------------------------
    std::vector<int> sendCounts(size, 0);

    for (int r = 0; r < nbinfoSize; ++r) {
        const int base = r * nbRecordSize;

        int eLoc       = nbinfo[base + 0];
        int lf_self    = nbinfo[base + 1]; (void)lf_self;
        int selfRank   = nbinfo[base + 2];
        int eglonb     = nbinfo[base + 3]; (void)eglonb;
        int lfnb       = nbinfo[base + 4]; (void)lfnb;
        int nghRank    = nbinfo[base + 5];

        if (selfRank != rank) {
            // nbinfo must have selfRank == this rank, but we ignore if not.
            continue;
        }
        if (nghRank < 0 || nghRank >= size) {
            continue;
        }

        int ifaceIdx = (eLoc >= 0 && eLoc < ne_local) ? local2iface[eLoc] : -1;
        if (ifaceIdx < 0) {
            // This local element is not in interfaceLocal (shouldn't happen if
            // nbinfo and classification are consistent).
            continue;
        }

        // One ifaceRecord per nbinfo record goes to nghRank
        sendCounts[nghRank] += ifaceRecordSize;
    }

    // ------------------------------------------------------------------
    // 3. Compute displacements and total send size, exchange counts
    // ------------------------------------------------------------------
    std::vector<int> sendDispls(size);
    int totalSendInts = 0;
    prefixSums(sendCounts, sendDispls, totalSendInts);

    recvCounts.assign(size, 0);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT,
                 recvCounts.data(), 1, MPI_INT,
                 comm);

    int totalRecvInts = 0;
    prefixSums(recvCounts, recvDispls, totalRecvInts);

    // ------------------------------------------------------------------
    // 4. Pack send buffer
    // ------------------------------------------------------------------
    std::vector<int> sendBuf(totalSendInts);
    std::vector<int> faceCounter(size, 0); // counts in "records" per dest

    for (int r = 0; r < nbinfoSize; ++r) {
        const int base = r * nbRecordSize;

        int eLoc       = nbinfo[base + 0];
        int lf_self    = nbinfo[base + 1];
        int selfRank   = nbinfo[base + 2];
        int eglonb     = nbinfo[base + 3]; // global elem ID on neighbor
        int lfnb       = nbinfo[base + 4]; // local face index on neighbor
        int nghRank    = nbinfo[base + 5];

        if (selfRank != rank) continue;
        if (nghRank < 0 || nghRank >= size) continue;

        int ifaceIdx = (eLoc >= 0 && eLoc < ne_local) ? local2iface[eLoc] : -1;
        if (ifaceIdx < 0) continue;

        int recordIdxInDest = faceCounter[nghRank]++; // 0-based record index
        int dstBase = sendDispls[nghRank] + recordIdxInDest * ifaceRecordSize;

        // Fill: [0] = eglonb, [1] = lfnb
        sendBuf[dstBase + 0] = eglonb;
        sendBuf[dstBase + 1] = lfnb;

        // Row index in ifaceNeigh arrays
        const int rowOffset = ifaceIdx * nfe;

        // [2 .. 2+nfe-1] = ifaceNeighGlobal row
        for (int lf = 0; lf < nfe; ++lf) {
            sendBuf[dstBase + 2 + lf] = ifaceNeighGlobal[rowOffset + lf];
        }
        // [2+nfe .. 2+2*nfe-1] = ifaceNeighOwner row
        for (int lf = 0; lf < nfe; ++lf) {
            sendBuf[dstBase + 2 + nfe + lf] = ifaceNeighOwner[rowOffset + lf];
        }
    }

    // ------------------------------------------------------------------
    // 5. Exchange via Alltoallv
    // ------------------------------------------------------------------
    recvBuf.resize(totalRecvInts);

    MPI_Alltoallv(sendBuf.data(),
                  sendCounts.data(), sendDispls.data(), MPI_INT,
                  recvBuf.data(),
                  recvCounts.data(), recvDispls.data(), MPI_INT,
                  comm);
}


void exchangeInterfaceNeighborRowsP2P(const Mesh&                  mesh,
                                   const ElementClassification& cls,
                                   const std::vector<int>&      ifaceNeighGlobal,
                                   const std::vector<int>&      ifaceNeighOwner,
                                   const int*                   nbinfo,
                                   int                          nbinfoSize,
                                   int                          nfe,
                                   MPI_Comm                     comm,
                                   std::vector<int>&            recvBuf,
                                   std::vector<int>&            recvCounts,
                                   std::vector<int>&            recvDispls)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    const int ne_local       = mesh.ne;
    const int nbRecordSize   = 6;            // nbinfo record size
    const int ifaceRecordSize = 2 + 2*nfe;   // eglonb, lfnb, G[nfe], Owner[nfe]

    // Quick exit if nothing to do
    if (ne_local == 0 || cls.interfaceLocal.empty() || nbinfoSize == 0) {
        recvBuf.clear();
        recvCounts.assign(size, 0);
        recvDispls.assign(size, 0);
        return;
    }

    // ------------------------------------------------------------------
    // 1. Build mapping: local element -> interface index
    // ------------------------------------------------------------------
    std::vector<int> local2iface(ne_local, -1);
    for (int k = 0; k < (int)cls.interfaceLocal.size(); ++k) {
        int eLoc = cls.interfaceLocal[k];
        if (eLoc >= 0 && eLoc < ne_local) {
            local2iface[eLoc] = k;
        }
    }

    // ------------------------------------------------------------------
    // 2. Count how many *ints* we will send to each rank, and collect
    //    the list of actual neighbor ranks (those with nonzero count).
    // ------------------------------------------------------------------
    std::vector<int> sendCounts(size, 0);  // counts in ints

    for (int r = 0; r < nbinfoSize; ++r) {
        const int base = r * nbRecordSize;

        int eLoc       = nbinfo[base + 0];
        int lf_self    = nbinfo[base + 1]; (void)lf_self;
        int selfRank   = nbinfo[base + 2];
        int eglonb     = nbinfo[base + 3]; (void)eglonb;
        int lfnb       = nbinfo[base + 4]; (void)lfnb;
        int nghRank    = nbinfo[base + 5];

        if (selfRank != rank) continue;
        if (nghRank < 0 || nghRank >= size) continue;

        int ifaceIdx = (eLoc >= 0 && eLoc < ne_local) ? local2iface[eLoc] : -1;
        if (ifaceIdx < 0) continue;

        // One ifaceRecord (ifaceRecordSize ints) per nbinfo record
        sendCounts[nghRank] += ifaceRecordSize;
    }

    // Build neighbor rank list = ranks with sendCounts[r] > 0
    std::vector<int> neighbors;
    neighbors.reserve(size);
    for (int r = 0; r < size; ++r) {
        if (sendCounts[r] > 0) {
            neighbors.push_back(r);
        }
    }

    if (neighbors.empty()) {
        // No neighbors to talk to
        recvBuf.clear();
        recvCounts.assign(size, 0);
        recvDispls.assign(size, 0);
        return;
    }

    // ------------------------------------------------------------------
    // 3. Compute displacements for send buffer, total send size
    // ------------------------------------------------------------------
    std::vector<int> sendDispls(size);
    int totalSendInts = 0;
    prefixSums(sendCounts, sendDispls, totalSendInts);

    // ------------------------------------------------------------------
    // 4. Pack send buffer (single contiguous buffer, but only parts for neighbors used)
    // ------------------------------------------------------------------
    std::vector<int> sendBuf(totalSendInts);
    std::vector<int> faceCounter(size, 0); // counts in records per dest

    for (int r = 0; r < nbinfoSize; ++r) {
        const int base = r * nbRecordSize;

        int eLoc       = nbinfo[base + 0];
        int lf_self    = nbinfo[base + 1];
        int selfRank   = nbinfo[base + 2];
        int eglonb     = nbinfo[base + 3]; // global elem ID on neighbor
        int lfnb       = nbinfo[base + 4]; // local face index on neighbor
        int nghRank    = nbinfo[base + 5];

        if (selfRank != rank) continue;
        if (nghRank < 0 || nghRank >= size) continue;

        int ifaceIdx = (eLoc >= 0 && eLoc < ne_local) ? local2iface[eLoc] : -1;
        if (ifaceIdx < 0) continue;

        int recordIdxInDest = faceCounter[nghRank]++;      // 0-based record index
        int dstBase         = sendDispls[nghRank] + recordIdxInDest * ifaceRecordSize;

        // [0] = eglonb, [1] = lfnb
        sendBuf[dstBase + 0] = eglonb;
        sendBuf[dstBase + 1] = lfnb;

        const int rowOffset = ifaceIdx * nfe;

        // [2 .. 2+nfe-1] = ifaceNeighGlobal row
        for (int lf = 0; lf < nfe; ++lf) {
            sendBuf[dstBase + 2 + lf] = ifaceNeighGlobal[rowOffset + lf];
        }
        // [2+nfe .. 2+2*nfe-1] = ifaceNeighOwner row
        for (int lf = 0; lf < nfe; ++lf) {
            sendBuf[dstBase + 2 + nfe + lf] = ifaceNeighOwner[rowOffset + lf];
        }
    }

    // ------------------------------------------------------------------
    // 5. Point-to-point exchange of counts (phase 1)
    // ------------------------------------------------------------------
    const int tagCounts = 9101;  // arbitrary but fixed tags
    const int tagData   = 9102;

    const int numNeighbors = static_cast<int>(neighbors.size());

    // For each neighbor we know how many ints we send → sendCounts[neighbor]
    // We need to receive counts from each neighbor.
    std::vector<int> recvCountsNeighbor(numNeighbors, 0);

    std::vector<MPI_Request> reqsCountSend(numNeighbors);
    std::vector<MPI_Request> reqsCountRecv(numNeighbors);

    for (int i = 0; i < numNeighbors; ++i) {
        int nbRank = neighbors[i];
        int sendCountInts = sendCounts[nbRank];

        // Send our count to nbRank
        MPI_Isend(&sendCountInts, 1, MPI_INT,
                  nbRank, tagCounts, comm, &reqsCountSend[i]);

        // Receive their count
        MPI_Irecv(&recvCountsNeighbor[i], 1, MPI_INT,
                  nbRank, tagCounts, comm, &reqsCountRecv[i]);
    }

    // Wait for all count exchanges
    if (numNeighbors > 0) {
        MPI_Waitall(numNeighbors, reqsCountSend.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(numNeighbors, reqsCountRecv.data(), MPI_STATUSES_IGNORE);
    }

    // ------------------------------------------------------------------
    // 6. Build recvCounts/recvDispls for *all* ranks (length = size),
    //    and compute total recv buffer size.
    // ------------------------------------------------------------------
    recvCounts.assign(size, 0);
    recvDispls.assign(size, 0);

    int totalRecvInts = 0;
    for (int i = 0; i < numNeighbors; ++i) {
        int nbRank = neighbors[i];
        recvCounts[nbRank] = recvCountsNeighbor[i];
    }

    // We want recvDispls to be consistent with recvCounts and recvBuf
    for (int r = 0; r < size; ++r) {
        recvDispls[r] = totalRecvInts;
        totalRecvInts += recvCounts[r];
    }

    recvBuf.resize(totalRecvInts);

    // ------------------------------------------------------------------
    // 7. Point-to-point exchange of actual data (phase 2)
    // ------------------------------------------------------------------
    std::vector<MPI_Request> reqsDataSend(numNeighbors);
    std::vector<MPI_Request> reqsDataRecv(numNeighbors);

    for (int i = 0; i < numNeighbors; ++i) {
        int nbRank = neighbors[i];

        int sendCountInts = sendCounts[nbRank];
        int recvCountInts = recvCounts[nbRank];

        int* sendPtr = (sendCountInts > 0)
                       ? sendBuf.data() + sendDispls[nbRank]
                       : nullptr;
        int* recvPtr = (recvCountInts > 0)
                       ? recvBuf.data() + recvDispls[nbRank]
                       : nullptr;

        // Post receive first
        if (recvCountInts > 0) {
            MPI_Irecv(recvPtr, recvCountInts, MPI_INT,
                      nbRank, tagData, comm, &reqsDataRecv[i]);
        } else {
            reqsDataRecv[i] = MPI_REQUEST_NULL;
        }

        // Then send
        if (sendCountInts > 0) {
            MPI_Isend(sendPtr, sendCountInts, MPI_INT,
                      nbRank, tagData, comm, &reqsDataSend[i]);
        } else {
            reqsDataSend[i] = MPI_REQUEST_NULL;
        }
    }

    if (numNeighbors > 0) {
        MPI_Waitall(numNeighbors, reqsDataRecv.data(), MPI_STATUSES_IGNORE);
        MPI_Waitall(numNeighbors, reqsDataSend.data(), MPI_STATUSES_IGNORE);
    }
}

// -----------------------------------------------------------------------------
// Step 4: on receiver side, use recvBuf to build outerElemGlobal/outerElemRank.
//
//   recvBuf records (per record):
//     [0]         = eglonb    (global elem ID on this rank)
//     [1]         = lfnb      (local face index on this rank)
//     [2 .. 2+nfe-1]         = neighbors' global IDs (of the remote element)
//     [2+nfe .. 2+2*nfe-1]   = neighbors' owner ranks
//
// We define outer elements as:
//
//   - neighbors of neighbor elements
//   - that do NOT belong to the current rank
//   - and are NOT already in neighborElemGlobal/neighborElemRank
//
// The function updates cls.outerElemGlobal / cls.outerElemRank, deduplicated.
// -----------------------------------------------------------------------------
void updateOuterElementsFromRecvBuf(ElementClassification&      cls,
                                    const Mesh&                mesh,
                                    const std::vector<int>&    recvBuf,
                                    const std::vector<int>&    recvCounts,
                                    const std::vector<int>&    recvDispls,
                                    int                        nfe,
                                    MPI_Comm                   comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    (void)size;

    cls.outerElemGlobal.clear();
    cls.outerElemRank.clear();

    if (recvBuf.empty())
        return;

    const int ifaceRecordSize = 2 + 2*nfe;

    // Set of local elements' global IDs
    std::unordered_set<int> myElems;
    myElems.reserve(mesh.ne * 2);
    for (int e = 0; e < mesh.ne; ++e) {
        myElems.insert(mesh.elemGlobalID[e]);
    }

    // Set of direct neighbor elements (globalID, rank)
    struct NeighborKey {
        int gid;
        int rank;
        bool operator==(const NeighborKey& o) const noexcept {
            return gid == o.gid && rank == o.rank;
        }
    };
    struct NeighborKeyHash {
        std::size_t operator()(const NeighborKey& k) const noexcept {
            std::size_t h1 = std::hash<int>{}(k.gid);
            std::size_t h2 = std::hash<int>{}(k.rank);
            return h1 ^ (h2 + 0x9e3779b9 + (h1<<6) + (h1>>2));
        }
    };

    std::unordered_set<NeighborKey, NeighborKeyHash> neighborSet;
    neighborSet.reserve(cls.neighborElemGlobal.size() * 2);
    for (std::size_t i = 0; i < cls.neighborElemGlobal.size(); ++i) {
        neighborSet.insert(NeighborKey{cls.neighborElemGlobal[i],
                                       cls.neighborElemRank[i]});
    }

    // Set of outer elements (globalID, rank) to dedupe
    std::unordered_set<NeighborKey, NeighborKeyHash> outerSet;

    // Walk through received records
    for (int src = 0; src < size; ++src) {
        int baseInts  = recvDispls[src];
        int countInts = recvCounts[src];
        if (countInts == 0) continue;

        int numRecords = (ifaceRecordSize > 0) ? (countInts / ifaceRecordSize) : 0;

        for (int rIdx = 0; rIdx < numRecords; ++rIdx) {
            int base = baseInts + rIdx * ifaceRecordSize;

            int eglonb = recvBuf[base + 0]; // our own global elem ID
            int lfnb   = recvBuf[base + 1]; (void)lfnb;

            const int* neighGlobalRow = &recvBuf[base + 2];
            const int* neighOwnerRow  = &recvBuf[base + 2 + nfe];

            // For each neighbor-of-neighbor
            for (int lf2 = 0; lf2 < nfe; ++lf2) {
                int candG = neighGlobalRow[lf2];
                int candR = neighOwnerRow[lf2];

                if (candG < 0)   continue;       // no element / boundary
                if (candR < 0)   continue;       // unknown owner
                if (candR == rank) continue;     // belongs to this rank → not "outer"
                if (candG == eglonb) continue;   // it's just us

                // skip if this is already a direct neighbor element
                NeighborKey nk{candG, candR};
                if (neighborSet.find(nk) != neighborSet.end())
                    continue;

                // skip if it's actually one of our local elements
                if (myElems.find(candG) != myElems.end())
                    continue;

                // True outer element
                if (outerSet.insert(nk).second) {
                    cls.outerElemGlobal.push_back(candG);
                    cls.outerElemRank.push_back(candR);
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// High-level driver: from e2e + nbinfo, update outerElemGlobal/outerElemRank
// in ElementClassification (assuming cls is already filled with interior,
// boundary, interface, neighbor elements).
// -----------------------------------------------------------------------------
void updateOuterElements(ElementClassification& cls,
                         const Mesh&           mesh,
                         const int*            e2e,        // [mesh.ne * nfe], global IDs or -1
                         int                   nfe,
                         const int*            nbinfo,     // 6*nbinfoSize ints
                         int                   nbinfoSize,
                         MPI_Comm              comm)
{
    // Step 1: gather neighbor globals for interface elements
    std::vector<int> ifaceNeighGlobal;
    collectInterfaceNeighborGlobals(e2e, nfe, cls, ifaceNeighGlobal);

    // Step 2: determine owner rank for each global ID (local + neighbor info only)
    std::vector<int> ifaceNeighOwner;
    determineOwnersForInterfaceNeighborsLocal(mesh, cls,
                                              ifaceNeighGlobal,
                                              ifaceNeighOwner,
                                              comm);

    // Step 3: exchange appropriate iface rows with neighbors
    std::vector<int> recvBuf, recvCounts, recvDispls;
    exchangeInterfaceNeighborRowsP2P(mesh, cls,
                                  ifaceNeighGlobal, ifaceNeighOwner,
                                  nbinfo, nbinfoSize, nfe,
                                  comm,
                                  recvBuf, recvCounts, recvDispls);

    // Step 4: on the receiver side, build outer elements
    updateOuterElementsFromRecvBuf(cls, mesh,
                                   recvBuf, recvCounts, recvDispls,
                                   nfe, comm);
}

DMD initializeDMD(Mesh& mesh, const Master& master, const PDE& pde, MPI_Comm comm) 
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  
    DMD dmd;

    mke2e_fill_first_neighbors(mesh.t2t.data(), mesh.t.data(), mesh.localfaces.data(),
               mesh.elemGlobalID.data(), mesh.nodeGlobalID.data(), mesh.ne, mesh.nve, 
              mesh.nvf, mesh.nfe, comm, dmd.nbinfo);
    
    int nboufaces = setboundaryfaces(mesh.t2t.data(), mesh.t.data(), mesh.localfaces.data(), mesh.p.data(),    
        mesh.boundaryExprs, mesh.dim, mesh.nve, mesh.nvf, mesh.nfe, mesh.ne, mesh.nbndexpr);

    if (pde.xdgfile == "") {
        mesh.xdg.resize(master.npe*mesh.dim*mesh.ne);
        compute_dgnodes(mesh.xdg.data(), mesh.p.data(), mesh.t.data(), master.phielem.data(), master.npe, mesh.dim, mesh.ne, mesh.nve);      
        project_dgnodes_onto_curved_boundaries(mesh.xdg.data(), mesh.t2t.data(), master.perm.data(), mesh.curvedBoundaries.data(),
                mesh.curvedBoundaryExprs, mesh.dim, master.porder, master.npe, master.npf, mesh.nfe, mesh.ne, -1);               
        if (rank==0) std::cout << "Finished computing dgnodes.\n";
    }     
    else 
      readParFieldFromBinaryFile(pde.xdgfile, mesh.elemGlobalID, mesh.xdg, mesh.xdgdims);
          
    setperiodicfaces(mesh.t2t.data(), mesh.t.data(), mesh.localfaces.data(), mesh.p.data(),    
        mesh.elemGlobalID.data(), mesh.periodicBoundaries1.data(), mesh.periodicBoundaries2.data(),
        mesh.periodicExprs1, mesh.periodicExprs2, mesh.dim, mesh.nve, mesh.nvf, mesh.nfe, mesh.ne, 
        mesh.nprdexpr, mesh.nprdcom, nboufaces, comm, dmd.nbinfo);
                                        
    dmd.numneigh = static_cast<int>(dmd.nbinfo.size() / 6); 
    ElementClassification elemclass;      
    classifyElementsWithE2EAndNbinfo(mesh.t2t.data(), mesh.elemGlobalID.data(), mesh.nfe, 
                                     mesh.ne, dmd.nbinfo.data(), dmd.numneigh, elemclass, rank);
    
    if (pde.hybrid==1) { // HDG
      elemclass.outerElemLocal.clear();
      elemclass.outerElemGlobal.clear();
      elemclass.outerElemRank.clear();
    }
    else if (pde.hybrid==0) { // LDG
      updateOuterElements(elemclass, mesh, mesh.t2t.data(), mesh.nfe, 
                            dmd.nbinfo.data(), dmd.numneigh,  comm);
    }
    
    buildElempartFromClassification(elemclass, dmd, rank);
    buildElem2CpuFromClassification(elemclass, dmd, comm);      
    buildElemRecv(dmd, comm);
    buildElemsend(mesh, dmd, comm);
          
    int nsend = dmd.elemsend.size();
    int nrecv = dmd.elemrecv.size();
    dmd.localelemsend.resize(nsend); 
    dmd.localelemrecv.resize(nrecv);
    for (int i=0; i<nsend; i++) dmd.localelemsend[i] = dmd.elemsend[i][1];
    for (int i=0; i<nrecv; i++) dmd.localelemrecv[i] = dmd.elemrecv[i][1];      

    return dmd;
}

void writemesh(Mesh& mesh, const DMD& dmd, const PDE& pde, const Master& master, MPI_Comm comm)
{    
    int nve = mesh.nve;
    int nfe = mesh.nfe;
    int ne = dmd.elempart.size();      
    mesh.bf.resize(nfe * ne);
    for (int i = 0; i < mesh.ne; i++) {
      int k = dmd.elempart_local[i];
      for (int j = 0; j < mesh.nfe; j++)
        mesh.bf[j + mesh.nfe*i] = (mesh.t2t[j + mesh.nfe*k] < 0) ? -mesh.t2t[j + mesh.nfe*k] : 0;
    }            
    // apply boundary conditions    
    for (int i = 0; i < mesh.nfe*mesh.ne; ++i) 
        if (mesh.bf[i] > 0) mesh.bf[i] = mesh.boundaryConditions[mesh.bf[i]-1];
        
    sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, mesh.bf, mesh.bf, nfe);

    computeElementToGlobalNodeMap(mesh, dmd.elempart_local);
    mesh.tg.resize(nve * ne);
    sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, mesh.tg, mesh.tg, mesh.nve);
  
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::string filename = make_path(pde.datainpath, "mesh" + std::to_string(rank+1) + ".bin");
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) error("Error opening file: " + filename);
    
    std::vector<int> ndims(20, 0);
    ndims[0] = (mesh.dim);
    ndims[1] = (dmd.elempart.size());
    ndims[3] = mesh.np_global;  
    ndims[4] = (mesh.nfe);
     
    std::vector<int> nsize(50, 0);
    nsize[0]  = 20; 
    if (size > 1) {
      nsize[4]  = (dmd.nbsd.size());
      nsize[5]  = (dmd.localelemsend.size());
      nsize[6]  = (dmd.localelemrecv.size());
      nsize[7]  = (dmd.elemsendpts.size());
      nsize[8]  = (dmd.elemrecvpts.size());
    }
    nsize[9]  = (dmd.elempart.size());
    nsize[10] = (dmd.elempartpts.size());    
    

    nsize[23] = (master.perm.size());
    nsize[24] = (mesh.bf.size());
    nsize[25] = (mesh.cartGridPart.size());
        
    nsize[26] = ne * nve;
    nsize[27] = mesh.boundaryConditions.size();
    nsize[28] = dmd.intepartpts.size();

    writeDouble(out, static_cast<double>(nsize.size())); 
    writeVectorAsDoubles(out, nsize);
    writeVectorAsDoubles(out, ndims);

    if (size > 1) {            
      writeVectorAsDoubles(out, dmd.nbsd);
      writeVectorAsDoubles(out, dmd.localelemsend);
      writeVectorAsDoubles(out, dmd.localelemrecv);
      writeVectorAsDoubles(out, dmd.elemsendpts);
      writeVectorAsDoubles(out, dmd.elemrecvpts);
    }
    writeVectorAsDoubles(out, dmd.elempart);
    writeVectorAsDoubles(out, dmd.elempartpts);

    if (pde.hybrid > 0) {
        writeVectorAsDoubles(out, master.perm);  
        writeVectorAsDoubles(out, mesh.bf);      
        writeVectorAsDoubles(out, mesh.cartGridPart);
    }

    writeVectorAsDoubles(out, mesh.tg);
    writeVectorAsDoubles(out, mesh.boundaryConditions);
    writeVectorAsDoubles(out, dmd.intepartpts);

    out.close();    
    std::cout << "Finished writing mesh to " + filename << std::endl;
}

void writesol(Mesh& mesh, const DMD& dmd, const PDE& pde, const Master& master, MPI_Comm comm)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
  
    std::string filename = make_path(pde.datainpath, "sol" + std::to_string(rank+1) + ".bin");
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) error("Error opening file: " + filename);
    
    int ne = dmd.elempart.size();     

    std::vector<double> ndims(12, 0);
    ndims[0] = ne;    
    ndims[2] = mesh.nfe;
    ndims[3] = master.npe;
    ndims[4] = master.npf;
    ndims[5] = pde.nc;
    ndims[6] = pde.ncu;
    ndims[7] = pde.ncq;
    ndims[8] = pde.ncw;
    ndims[9] = pde.ncv;
    ndims[10] = pde.nch;
    ndims[11] = pde.ncx;    

    std::vector<double> nsize(20, 0);
    nsize[0] = static_cast<double>(ndims.size());
    nsize[1] = master.npe*mesh.dim*ne;

    if (pde.udgfile != "") {
      readParFieldFromBinaryFile(pde.udgfile, mesh.elemGlobalID, mesh.udg, mesh.udgdims);      
      nsize[2] = master.npe*mesh.udgdims[1]*ne;
    }
    if (pde.vdgfile != "") {
      readParFieldFromBinaryFile(pde.vdgfile, mesh.elemGlobalID, mesh.vdg, mesh.vdgdims);   
      nsize[3] = master.npe*mesh.vdgdims[1]*ne;
    }
    if (pde.wdgfile != "") {
      readParFieldFromBinaryFile(pde.wdgfile, mesh.elemGlobalID, mesh.wdg, mesh.wdgdims);   
      nsize[4] = master.npe*mesh.wdgdims[1]*ne;
    }

    // Helper to write vector as binary double
    auto writeDoubleVector = [&](const std::vector<double>& vec) {
        out.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));
    };

    double len = static_cast<double>(nsize.size());
    out.write(reinterpret_cast<const char*>(&len), sizeof(double));
    writeDoubleVector(nsize);
    writeDoubleVector(ndims);

    vector<double> xdg(master.npe*mesh.dim*ne, 0);
    select_columns(xdg.data(), mesh.xdg.data(), dmd.elempart_local.data(), master.npe*mesh.dim, mesh.ne);
    sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*mesh.dim);
    writeDoubleVector(xdg);
    
    if (pde.udgfile != "") {
      int nc = mesh.udgdims[1];
      xdg.resize(master.npe*nc*ne, 0);
      select_columns(xdg.data(), mesh.udg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      writeDoubleVector(xdg);
    }
    if (pde.vdgfile != "") {
      int nc = mesh.vdgdims[1];
      xdg.resize(master.npe*nc*ne, 0);
      select_columns(xdg.data(), mesh.vdg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      writeDoubleVector(xdg);
    }
    if (pde.wdgfile != "") {
      int nc = mesh.wdgdims[1];
      xdg.resize(master.npe*nc*ne, 0);
      select_columns(xdg.data(), mesh.wdg.data(), dmd.elempart_local.data(), master.npe*nc, mesh.ne);
      sendrecvdata(comm, dmd.nbsd, dmd.elemsendpts, dmd.elemrecvpts, 
                 dmd.localelemsend, dmd.localelemrecv, xdg, xdg, master.npe*nc);
      writeDoubleVector(xdg);
    }
      
    out.close();    
    std::cout << "Finished writing initial solution to " + filename << std::endl;
}

// // nbinfo: output array of (6 + nfe) ints per remote face:
// //   [0]            = eLoc        (local element index on this rank)
// //   [1]            = lf          (local face index on this rank)
// //   [2]            = self_rank   (this rank)
// //   [3]            = eglonb      (neighbor global element ID)
// //   [4]            = lfnb        (neighbor local face index on neighbor_rank)
// //   [5]            = neighbor_rank
// //   [6..6+nfe-1]   = e2e_nb[:]   (neighbor's e2e row for that element, global IDs)
// // nbinfoSize: number of such records (so total ints = nbinfoSize * (6 + nfe)).
// //
// // e2e on entry is assumed to already contain local neighbors (global IDs)
// // and -1 for remote/boundary faces. This function fills remote faces and nbinfo.
// void mke2e_fill_second_neighbors(int*       e2e,
//                                  const int* e2n,
//                                  const int* local_faces,
//                                  const int* elemGlobalID,
//                                  const int* nodeGlobalID,
//                                  int        ne,
//                                  int        nne,
//                                  int        nnf,
//                                  int        nfe,
//                                  MPI_Comm   comm,
//                                  int*&      nbinfo,
//                                  int&       nbinfoSize)
// {
//     int rank, size;
//     MPI_Comm_rank(comm, &rank);
//     MPI_Comm_size(comm, &size);
// 
//     nbinfo     = nullptr;
//     nbinfoSize = 0;
// 
//     if (nnf > 4) {
//         error("mke2e_fill_remote_neighbors: nnf > 4 not supported in this implementation");
//     }
// 
//     // ------------------------------------------------------------------
//     // 1. Collect all faces with e2e == -1 (potential boundary or remote)
//     //    Owner rule: owner = min(global_node_ids) % size
//     // ------------------------------------------------------------------
//     struct LocalFaceRecord {
//         std::array<int,4> gnodes; // global node IDs of this face
//         int               elemGlobal;
//         int               elemLocal;
//         int               lf;
//     };
// 
//     std::vector<int>             sendFaceCounts(size, 0);
//     std::vector<LocalFaceRecord> exportFaces;
//     exportFaces.reserve(ne * nfe);
// 
//     for (int e = 0; e < ne; ++e) {
//         int eGlob = elemGlobalID[e];
//         for (int lf = 0; lf < nfe; ++lf) {
//             if (e2e[e * nfe + lf] != -1) {
//                 // already has on-rank neighbor, skip
//                 continue;
//             }
// 
//             LocalFaceRecord rec;
//             rec.elemGlobal = eGlob;
//             rec.elemLocal  = e;
//             rec.lf         = lf;
// 
//             int minNode = std::numeric_limits<int>::max();
//             for (int i = 0; i < nnf; ++i) {
//                 int ln  = local_faces[lf * nnf + i];  // local node in element
//                 int ln2 = e2n[e * nne + ln];          // local node in mesh
//                 int gn  = nodeGlobalID[ln2];          // global node ID
//                 rec.gnodes[i] = gn;
//                 if (gn < minNode) minNode = gn;
//             }
// 
//             int owner = (minNode % size + size) % size;
//             sendFaceCounts[owner] += 1;
//             exportFaces.push_back(rec);
//         }
//     }
// 
//     // ------------------------------------------------------------------
//     // 2. Pack faces into a send buffer and Alltoallv to owners
//     //
//     // For each face record we also send the full e2e row of that element.
//     //
//     // Layout (ints):
//     //   [0..nnf-1]     : global node IDs
//     //   [nnf]          : elemGlobal
//     //   [nnf+1]        : elemLocal
//     //   [nnf+2]        : lf
//     //   [nnf+3 .. nnf+3+nfe-1] : e2eRow (global neighbor IDs of this element)
//     //
//     // recordSize = nnf + 3 + nfe
//     // ------------------------------------------------------------------
//     const int recordSize = nnf + 3 + nfe;
// 
//     std::vector<int> sendCountsInts(size), sendDisplsInts;
//     std::vector<int> recvCountsInts(size), recvDisplsInts;
//     for (int r = 0; r < size; ++r) {
//         sendCountsInts[r] = sendFaceCounts[r] * recordSize;
//     }
// 
//     int totalSendInts = 0, totalRecvInts = 0;
//     prefixSums(sendCountsInts, sendDisplsInts, totalSendInts);
// 
//     MPI_Alltoall(sendCountsInts.data(), 1, MPI_INT,
//                  recvCountsInts.data(), 1, MPI_INT,
//                  comm);
// 
//     prefixSums(recvCountsInts, recvDisplsInts, totalRecvInts);
// 
//     std::vector<int> sendBuf(totalSendInts);
//     std::vector<int> faceOffsetPerDest(size, 0); // in faces
// 
//     for (const auto& rec : exportFaces) {
//         int minNode = std::numeric_limits<int>::max();
//         for (int i = 0; i < nnf; ++i)
//             if (rec.gnodes[i] < minNode) minNode = rec.gnodes[i];
// 
//         int owner = (minNode % size + size) % size;
// 
//         int faceIdxInDest = faceOffsetPerDest[owner]++;
//         int base = sendDisplsInts[owner] + faceIdxInDest * recordSize;
// 
//         // nodes
//         for (int i = 0; i < nnf; ++i) {
//             sendBuf[base + i] = rec.gnodes[i];
//         }
//         // element info
//         sendBuf[base + nnf    ] = rec.elemGlobal;
//         sendBuf[base + nnf + 1] = rec.elemLocal;
//         sendBuf[base + nnf + 2] = rec.lf;
//         // e2e row (global neighbors) for this local element
//         int eLoc = rec.elemLocal;
//         for (int f = 0; f < nfe; ++f) {
//             sendBuf[base + nnf + 3 + f] = e2e[eLoc * nfe + f];
//         }
//     }
// 
//     std::vector<int> recvBuf(totalRecvInts);
// 
//     MPI_Alltoallv(sendBuf.data(),
//                   sendCountsInts.data(), sendDisplsInts.data(), MPI_INT,
//                   recvBuf.data(),
//                   recvCountsInts.data(), recvDisplsInts.data(), MPI_INT,
//                   comm);
// 
//     // ------------------------------------------------------------------
//     // 3. On owner ranks: match faces by FaceKey and build neighbor info
//     // ------------------------------------------------------------------
//     struct FaceOwnerRec {
//         int              rank;        // original rank
//         int              elemGlobal;  // global element ID
//         int              elemLocal;   // local element on that rank
//         int              lf;          // local face index on that element
//         std::vector<int> e2eRow;      // neighbor row for this element (size nfe)
//     };
// 
//     std::unordered_map<FaceKey, std::vector<FaceOwnerRec>, FaceKeyHash> faceMap;
//     faceMap.reserve(static_cast<std::size_t>(totalRecvInts / recordSize) * 2);
// 
//     // Decode recvBuf into faceMap
//     for (int src = 0; src < size; ++src) {
//         int countInts = recvCountsInts[src];
//         int baseInts  = recvDisplsInts[src];
// 
//         int numFaces = (recordSize > 0) ? (countInts / recordSize) : 0;
//         for (int f = 0; f < numFaces; ++f) {
//             int base = baseInts + f * recordSize;
// 
//             FaceOwnerRec rec;
//             rec.rank       = src;
//             rec.elemGlobal = recvBuf[base + nnf];
//             rec.elemLocal  = recvBuf[base + nnf + 1];
//             rec.lf         = recvBuf[base + nnf + 2];
// 
//             rec.e2eRow.resize(nfe);
//             for (int ff = 0; ff < nfe; ++ff) {
//                 rec.e2eRow[ff] = recvBuf[base + nnf + 3 + ff];
//             }
// 
//             FaceKey key;
//             key.len = static_cast<uint8_t>(nnf);
//             for (int i = 0; i < nnf; ++i) {
//                 key.nodes[i] = recvBuf[base + i];
//             }
//             std::sort(key.nodes.begin(), key.nodes.begin() + nnf);
// 
//             auto& vec = faceMap[key];
//             vec.push_back(rec);
//         }
//     }
// 
//     // Updates to send back:
//     //   selfElemLocal, selfLf, neighGlobal, neighLf, neighRank, neighE2E[:]
//     struct Update {
//         int              selfElemLocal;
//         int              selfLf;
//         int              neighGlobal;
//         int              neighLf;
//         int              neighRank;
//         std::vector<int> neighE2E;    // size nfe
//     };
// 
//     std::vector<std::vector<Update>> updates(size);
// 
//     for (auto& kv : faceMap) {
//         auto& owners = kv.second;
//         if (owners.size() == 2) {
//             const auto& a = owners[0];
//             const auto& b = owners[1];
// 
//             // To rank a.rank: self = a, neighbor = b
//             {
//                 Update ua;
//                 ua.selfElemLocal = a.elemLocal;
//                 ua.selfLf        = a.lf;
//                 ua.neighGlobal   = b.elemGlobal;
//                 ua.neighLf       = b.lf;
//                 ua.neighRank     = b.rank;
//                 ua.neighE2E      = b.e2eRow;
//                 updates[a.rank].push_back(std::move(ua));
//             }
// 
//             // To rank b.rank: self = b, neighbor = a
//             {
//                 Update ub;
//                 ub.selfElemLocal = b.elemLocal;
//                 ub.selfLf        = b.lf;
//                 ub.neighGlobal   = a.elemGlobal;
//                 ub.neighLf       = a.lf;
//                 ub.neighRank     = a.rank;
//                 ub.neighE2E      = a.e2eRow;
//                 updates[b.rank].push_back(std::move(ub));
//             }
//         }
//         // owners.size() == 1 → boundary, nothing to do
//     }
// 
//     // ------------------------------------------------------------------
//     // 4. Send updates back, fill e2e and nbinfo
//     //
//     // Each update record (ints):
//     //   [0]                = selfElemLocal
//     //   [1]                = selfLf
//     //   [2]                = neighGlobal
//     //   [3]                = neighLf
//     //   [4]                = neighRank
//     //   [5..5+nfe-1]       = neighE2E[:]
//     //
//     // updRecordSize = 5 + nfe
//     // ------------------------------------------------------------------
//     const int updRecordSize = 5 + nfe;
// 
//     std::vector<int> updSendCountsInts(size), updSendDisplsInts;
//     std::vector<int> updRecvCountsInts(size), updRecvDisplsInts;
// 
//     for (int r = 0; r < size; ++r) {
//         updSendCountsInts[r] = static_cast<int>(updates[r].size()) * updRecordSize;
//     }
// 
//     int totalUpdSendInts = 0, totalUpdRecvInts = 0;
//     prefixSums(updSendCountsInts, updSendDisplsInts, totalUpdSendInts);
// 
//     MPI_Alltoall(updSendCountsInts.data(), 1, MPI_INT,
//                  updRecvCountsInts.data(), 1, MPI_INT,
//                  comm);
// 
//     prefixSums(updRecvCountsInts, updRecvDisplsInts, totalUpdRecvInts);
// 
//     std::vector<int> updSendBuf(totalUpdSendInts);
//     std::vector<int> updRecvBuf(totalUpdRecvInts);
// 
//     // Pack updates
//     for (int r = 0; r < size; ++r) {
//         int base = updSendDisplsInts[r];
//         const auto& ur = updates[r];
//         for (std::size_t k = 0; k < ur.size(); ++k) {
//             const Update& U = ur[k];
//             int b = base + static_cast<int>(k) * updRecordSize;
// 
//             updSendBuf[b + 0] = U.selfElemLocal;
//             updSendBuf[b + 1] = U.selfLf;
//             updSendBuf[b + 2] = U.neighGlobal;
//             updSendBuf[b + 3] = U.neighLf;
//             updSendBuf[b + 4] = U.neighRank;
//             for (int ff = 0; ff < nfe; ++ff) {
//                 updSendBuf[b + 5 + ff] = U.neighE2E[ff];
//             }
//         }
//     }
// 
//     MPI_Alltoallv(updSendBuf.data(),
//                   updSendCountsInts.data(), updSendDisplsInts.data(), MPI_INT,
//                   updRecvBuf.data(),
//                   updRecvCountsInts.data(), updRecvDisplsInts.data(), MPI_INT,
//                   comm);
// 
//     // Collect neighbor info: 6 + nfe ints per record
//     std::vector<int> nbinfoVec;
// 
//     const int nbRecordSize = 6 + nfe;
// 
//     // Apply updates locally and build nbinfo
//     for (int src = 0; src < size; ++src) {
//         int countInts = updRecvCountsInts[src];
//         int baseInts  = updRecvDisplsInts[src];
//         int numUpd    = (updRecordSize > 0) ? (countInts / updRecordSize) : 0;
// 
//         for (int k = 0; k < numUpd; ++k) {
//             int base      = baseInts + k * updRecordSize;
//             int eLoc      = updRecvBuf[base + 0];
//             int lf        = updRecvBuf[base + 1];
//             int nghGlobal = updRecvBuf[base + 2];
//             int lfnb      = updRecvBuf[base + 3];
//             int nghRank   = updRecvBuf[base + 4];
// 
//             // neighbor's e2e row
//             std::vector<int> neighRow(nfe);
//             for (int ff = 0; ff < nfe; ++ff) {
//                 neighRow[ff] = updRecvBuf[base + 5 + ff];
//             }
// 
//             // overwrite -1 with neighbor global element ID (for this face)
//             e2e[eLoc * nfe + lf] = nghGlobal;
// 
//             // Append one nbinfo record: (eLoc, lf, self_rank, eglonb, lfnb, neighbor_rank, e2e_nb[:])
//             nbinfoVec.reserve(nbinfoVec.size() + nbRecordSize);
// 
//             nbinfoVec.push_back(eLoc);
//             nbinfoVec.push_back(lf);
//             nbinfoVec.push_back(rank);      // self_rank
//             nbinfoVec.push_back(nghGlobal); // eglonb
//             nbinfoVec.push_back(lfnb);      // lfnb
//             nbinfoVec.push_back(nghRank);   // neighbor_rank
// 
//             for (int ff = 0; ff < nfe; ++ff) {
//                 nbinfoVec.push_back(neighRow[ff]); // e2e(:, elemGlobal_on_neighbor)
//             }
//         }
//     }
// 
//     // ------------------------------------------------------------------
//     // 5. Export nbinfo as raw array
//     // ------------------------------------------------------------------
//     if (!nbinfoVec.empty()) {
//         nbinfoSize = static_cast<int>(nbinfoVec.size() / nbRecordSize);
//         nbinfo = static_cast<int*>(std::malloc(nbinfoVec.size() * sizeof(int)));
//         if (!nbinfo) {
//             error("mke2e_fill_remote_neighbors: malloc failed for nbinfo");
//         }
//         std::memcpy(nbinfo, nbinfoVec.data(),
//                     nbinfoVec.size() * sizeof(int));
//     } else {
//         nbinfo     = nullptr;
//         nbinfoSize = 0;
//     }
// }

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

#endif