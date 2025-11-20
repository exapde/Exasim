class CVisualization {
public:
    float* scafields=nullptr;
    float* vecfields=nullptr;
    float* tenfields=nullptr;
    float* srffields=nullptr;

    // Geometry (owned; always 3D padded: [3 x npoints], column-major)
    std::vector<float>  cgnodes;
    
    // Topology (owned)
    std::vector<int32_t> cgcells;      // [nve x ncells], column-major (0-based)
    std::vector<int32_t> celloffsets;  // [ncells] cumulative ends (nve, 2*nve, ...)
    std::vector<uint8_t> celltypes;    // [ncells] VTK cell type per cell (uniform)

    // Geometry (owned; always 3D padded: [3 x nnodes], column-major)
    std::vector<float>  facenodes;
    
    // Topology (owned)
    std::vector<int32_t> faceconn;     // [nvf x nfaces], column-major (0-based)
    std::vector<int32_t> faceoffsets;  // [ncells] cumulative ends (nve, 2*nve, ...)
    std::vector<uint8_t> facetypes;    // [ncells] VTK cell type per cell (uniform)
    
    // Field schemas (names only; data are provided at write time)
    std::vector<std::string> scalar_names;   // nscalars
    std::vector<std::string> vector_names;   // nvectors (3 comps each, z padded if nd==2)
    std::vector<std::string> tensor_names;   // ntensors (ntc comps each, ntc=nd*nd)
    std::vector<std::string> surface_names;   // nsurfaces

    // how fields were allocated: 0=CPU malloc/free, 2=CUDA host (cudaHostAlloc),
    // 3=HIP  host (hipHostMalloc), anything else => unknown/none
    int host_alloc_backend = 0;

    // Sizes / meta
    int nd      = 0;   // spatial dimension (2 or 3)
    int ntc     = 0;   // tensor components = nd*nd (4 or 9)
    int npoints = 0;
    int ncells  = 0;
    int nve     = 0;
    int nnodes  = 0;
    int nfaces  = 0;
    int nvf     = 0;    
    int savemode = 0;
    int rank = 0;
    
    // Precomputed appended-data offsets for VTU metadata
    std::vector<std::uint64_t> scalar_offsets; // size = nscalars
    std::vector<std::uint64_t> vector_offsets; // size = nvectors
    std::vector<std::uint64_t> tensor_offsets; // size = ntensors
    std::uint64_t points_offset = 0;
    std::uint64_t conn_offset   = 0;
    std::uint64_t offs_offset   = 0;
    std::uint64_t types_offset  = 0;

public:
    CVisualization(const dstype* xcg, int nd_in, int np,
                   const int* cgelcon, int npe, int ne,
                   const int* telem,   int nce, int nverts_per_cell,
                   int elemtype,
                   const std::vector<std::string>& scalars,
                   const std::vector<std::string>& vectors,
                   const std::vector<std::string>& tensors,
                   const std::vector<std::string>& surfaces)
    {
        if (np > 0) {
            Init(xcg, nd_in, np, cgelcon, npe, ne,
                 telem, nce, nverts_per_cell, elemtype,
                 scalars, vectors, tensors, surfaces);            
        }
    }

    CVisualization(CDiscretization& disc, int backend) {      
        int rank = disc.common.mpiRank;
        int nd_in   = disc.common.nd;
        int npoints_in = disc.sol.szxcg / nd_in;
        
        if (npoints_in > 0) {            
            int porder  = disc.common.porder;        
            int nsca    = disc.common.nsca;
            int nvec    = disc.common.nvec;            
            int nten    = disc.common.nten;            
            int nsurf   = disc.common.nsurf;            
            int npe     = disc.common.npe;
            int ne      = disc.common.ne1;
            int elemtype= disc.common.elemtype;
            int nve_in  = (elemtype==0) ? (nd_in + 1) : std::pow(2, nd_in);
                
            std::string fn1 = make_path(disc.common.exasimpath, "text2code/text2code/masternodes.bin");
            std::vector<dstype> xpe, xpf;
            std::vector<int> telem, tface, perm;
            masternodes(xpe, telem, xpf, tface, perm, porder, nd_in, elemtype, fn1);
    
            int nce = (int)telem.size() / nve_in;
    
            std::vector<std::string> scalars(nsca);
            for (int i = 0; i < nsca; i++) scalars[i] = "Scalar Field " + std::to_string(i);
    
            std::vector<std::string> vectors(nvec);
            for (int i = 0; i < nvec; i++) vectors[i] = "Vector Field " + std::to_string(i);
    
            std::vector<std::string> tensors(nten);
            for (int i = 0; i < nten; i++) tensors[i] = "Tensor Field " + std::to_string(i);

            std::vector<std::string> surfaces(nsurf);
            for (int i = 0; i < nsurf; i++) surfaces[i] = "Surface Field " + std::to_string(i);
            
            int* cgelcon;
            if (backend==0) cgelcon = &disc.mesh.cgelcon[0];
            else {
                TemplateMalloc(&cgelcon, npe*ne, 0);
                TemplateCopytoHost(cgelcon, disc.mesh.cgelcon, npe*ne, backend);
            }            
            
            Init(disc.sol.xcg, nd_in, npoints_in, cgelcon, npe, ne,
                 telem.data(), nce, nve_in, elemtype,
                 scalars, vectors, tensors, surfaces);

            if (backend != 0) CPUFREE(cgelcon);    

            savemode = (nsca + nvec + nten > 0); 
        
            if (backend==2) { // GPU
            #ifdef HAVE_CUDA        
                cudaTemplateHostAlloc(&scafields, npoints*nsca, cudaHostAllocMapped); // zero copy
                cudaTemplateHostAlloc(&vecfields, 3*npoints*nvec, cudaHostAllocMapped); // zero copy
                cudaTemplateHostAlloc(&tenfields, ntc*npoints*nten, cudaHostAllocMapped); // zero copy
                host_alloc_backend = 2;
            #endif                  
            }
            else if (backend==3) { // GPU
            #ifdef HAVE_HIP        
                hipTemplateHostMalloc(&scafields, npoints*nsca, hipHostMallocMapped); // zero copy
                hipTemplateHostMalloc(&vecfields, 3*npoints*nvec, hipHostMallocMapped); // zero copy
                hipTemplateHostMalloc(&tenfields, ntc*npoints*nten, hipHostMallocMapped); // zero copy                
                host_alloc_backend = 3;
            #endif                  
            }    
            else { // CPU
                scafields = (float *) malloc(npoints*nsca*sizeof(float));
                vecfields = (float *) malloc(3*npoints*nvec*sizeof(float));
                tenfields = (float *) malloc(ntc*npoints*nten*sizeof(float));
                host_alloc_backend = 0;
            }
            
            //cout<<ne<<"  "<<npoints<<endl;
            if (disc.common.mpiRank == 0) printf("finish CVisualization constructor... \n");    
        }        
    }

    ~CVisualization() {
        // Free scafields / vecfields / tenfields / srffields according to how they were allocated.
        auto free_field = [this](float*& p) {
            if (!p) return;
            switch (host_alloc_backend) {
                case 2:  // CUDA pinned host
                #ifdef HAVE_CUDA
                    cudaFreeHost(p);
                    p = nullptr;
                    break;
                #endif
    
                case 3:  // HIP pinned host
                #ifdef HAVE_HIP
                    hipHostFree(p);
                    p = nullptr;
                    break;
                #endif
    
                case 0:  // CPU malloc
                     std::free(p);
                     p = nullptr;
                     break;

                default:
                    std::free(p);
                    p = nullptr;
                    break;
            }
        };

        free_field(scafields);
        free_field(vecfields);
        free_field(tenfields);
        free_field(srffields); // in case you allocate this later
        if (rank==0) printf("CVisualization is freed successfully.\n");
    }

    // ============================ Writers ============================

    // scalarfields: [npoints * nscalars]
    // vectorfields: [3*npoints * nvectors] (z padded if nd==2)
    // tensorfields: [ntc*npoints * ntensors], ntc = nd*nd (4 for 2D, 9 for 3D)
    void vtuwrite(const std::string& filename_no_ext,
                  const float* scalarfields,   // nullptr allowed iff scalar_names empty
                  const float* vectorfields = nullptr,   // nullptr allowed iff vector_names empty
                  const float* tensorfields = nullptr) const // nullptr allowed iff tensor_names empty
    {
        if (!scalar_names.empty() && !scalarfields)
            throw std::invalid_argument("vtuwrite: scalarfields is null but scalar_names is non-empty.");
        if (!vector_names.empty() && !vectorfields)
            throw std::invalid_argument("vtuwrite: vectorfields is null but vector_names is non-empty.");
        if (!tensor_names.empty() && !tensorfields)
            throw std::invalid_argument("vtuwrite: tensorfields is null but tensor_names is non-empty.");

        const std::string filename = filename_no_ext + ".vtu";
        std::ofstream os(filename, std::ios::binary);
        if (!os) throw std::runtime_error("Cannot open output file: " + filename);

        os << "<?xml version=\"1.0\"?>\n";
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\""
           << vtk_byte_order() << "\" header_type=\"UInt64\">\n";
        os << "  <UnstructuredGrid>\n";
        os << "    <Piece NumberOfPoints=\"" << npoints
           << "\" NumberOfCells=\"" << ncells << "\">\n";

        // PointData schema
        if (!scalar_names.empty() || !vector_names.empty() || !tensor_names.empty()) {
            os << "      <PointData Scalars=\"scalars\">\n";
            for (int i = 0; i < (int)scalar_names.size(); ++i)
                os << "        <DataArray type=\"Float32\" Name=\"" << scalar_names[i]
                   << "\" Format=\"appended\" offset=\"" << scalar_offsets[i] << "\"/>\n";
            for (int i = 0; i < (int)vector_names.size(); ++i)
                os << "        <DataArray type=\"Float32\" Name=\"" << vector_names[i]
                   << "\" NumberOfComponents=\"3\" Format=\"appended\" offset=\"" << vector_offsets[i] << "\"/>\n";
            for (int i = 0; i < (int)tensor_names.size(); ++i)
                os << "        <DataArray type=\"Float32\" Name=\"" << tensor_names[i]
                   << "\" NumberOfComponents=\"" << ntc << "\" Format=\"appended\" offset=\"" << tensor_offsets[i] << "\"/>\n";
            os << "      </PointData>\n";
        }

        // Points schema
        os << "      <Points>\n";
        os << "        <DataArray type=\"Float32\" Name=\"points\" NumberOfComponents=\"3\" Format=\"appended\" offset=\""
           << points_offset << "\"/>\n";
        os << "      </Points>\n";

        // Cells schema
        os << "      <Cells>\n";
        os << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"appended\" offset=\"" << conn_offset << "\"/>\n";
        os << "        <DataArray type=\"Int32\" Name=\"offsets\"     Format=\"appended\" offset=\"" << offs_offset << "\"/>\n";
        os << "        <DataArray type=\"UInt8\" Name=\"types\"       Format=\"appended\" offset=\"" << types_offset << "\"/>\n";
        os << "      </Cells>\n";

        os << "    </Piece>\n";
        os << "  </UnstructuredGrid>\n";
        os << "  <AppendedData encoding=\"raw\">\n";
        os << "   _";

        // Appended blocks
        for (int si = 0; si < (int)scalar_names.size(); ++si) {
            //const dstype* src = &scalarfields[npoints * si];
            //write_as_float32(os, src, npoints);
            write_block(os, &scalarfields[npoints * si], npoints * (int)sizeof(float));
        }       

        for (int vi = 0; vi < (int)vector_names.size(); ++vi) {
            //const dstype* src = &vectorfields[3 * npoints * vi];
            //write_as_float32(os, src, 3 * npoints);
            write_block(os, &vectorfields[3 * npoints * vi], 3 * npoints * (int)sizeof(float));
        }

        for (int ti = 0; ti < (int)tensor_names.size(); ++ti) {
            //const dstype* src = &tensorfields[ntc * npoints * ti];
            //write_as_float32(os, src, ntc * npoints);
            write_block(os, &tensorfields[ntc * npoints * ti], ntc * npoints * (int)sizeof(float));
        }

        write_block(os, cgnodes.data(),      3 * npoints * (int)sizeof(float));
        write_block(os, cgcells.data(),      nve * ncells * (int)sizeof(int32_t));
        write_block(os, celloffsets.data(),  ncells      * (int)sizeof(int32_t));
        write_block(os, celltypes.data(),    ncells      * (int)sizeof(uint8_t));

        os << "\n  </AppendedData>\n";
        os << "</VTKFile>\n";
        os.close();
    }

    // Parallel writer (rank pieces + PVTU on rank 0)
    void vtuwrite_parallel(const std::string& base_name,
                           int rank, int nranks,
                           const float* scalarfields,
                           const float* vectorfields = nullptr,
                           const float* tensorfields = nullptr) const
    {
        const std::string my_base = base_name + rank_tag(rank);
        vtuwrite(my_base, scalarfields, vectorfields, tensorfields);

        if (rank == 0) {
            std::vector<std::string> pieces;
            pieces.reserve(nranks);
            for (int r = 0; r < nranks; ++r)
                pieces.push_back(base_name + rank_tag(r) + ".vtu");
            write_pvtu(base_name, pieces, scalar_names, vector_names, tensor_names, ntc);
        }
    }

    // PVD writers (unchanged)
    static void pvdwrite(const std::string& pvd_name_no_ext,
                         const std::vector<std::string>& files,
                         const std::vector<float>& times) {
        if (files.size() != times.size())
            throw std::invalid_argument("pvdwrite: files and times must have same length.");
        const std::string pvdfile = pvd_name_no_ext + ".pvd";
        std::ofstream os(pvdfile, std::ios::binary);
        if (!os) throw std::runtime_error("Cannot open PVD file: " + pvdfile);
        os << "<?xml version=\"1.0\"?>\n";
        os << "<VTKFile type=\"Collection\" version=\"0.1\"\n";
        os << "         byte_order=\"" << vtk_byte_order() << "\"\n";
        os << "         compressor=\"vtkZLibDataCompressor\">\n";
        os << "  <Collection>\n";
        for (int i = 0; i < (int)files.size(); ++i) {
            os << "    <DataSet timestep=\"" << times[i] << "\" group=\"\" part=\"0\"\n";
            os << "             file=\"" << files[i] << "\"/>\n";
        }
        os << "  </Collection>\n";
        os << "</VTKFile>\n";
        if (!os) throw std::runtime_error("Error writing PVD file.");
    }

    static void pvdwrite_series(const std::string& base,
                                const dstype* dt, int nt, int nm,
                                const std::string& ext /* "vtu" or "pvtu" */) {
        if (!dt && nt > 0) throw std::invalid_argument("pvdwrite_series: dt pointer is null.");
        if (nt < 0)        throw std::invalid_argument("pvdwrite_series: nt must be non-negative.");
        if (!(ext == "vtu" || ext == "pvtu"))
            throw std::invalid_argument("pvdwrite_series: ext must be \"vtu\" or \"pvtu\".");
        std::vector<std::string> files; files.reserve(nt);
        std::vector<float>       times; times.reserve(nt);
        float t = 0.0f;
        for (int k = 0; k < nt; ++k) {
            t += dt[k];
            if ((k+1) % nm == 0) {
                files.push_back(base + "_" + step_tag(k+1) + "." + ext);
                times.push_back(t);
            }
        }
        cout<<nt<<",  "<<nm<<",  "<<base<<endl;        
        pvdwrite(base, files, times);
    }

private:
    void Init(const dstype* xcg, int nd_in, int np,
              const int* cgelcon, int npe, int ne,
              const int* telem,   int nce, int nverts_per_cell,
              int elemtype,
              const std::vector<std::string>& s_names,
              const std::vector<std::string>& v_names,
              const std::vector<std::string>& t_names,
              const std::vector<std::string>& surf_names)
    {
        nd = nd_in;
        ntc = nd * nd;           // <-- your requested change
        npoints = np;
        nve = nverts_per_cell;
        scalar_names = s_names;
        vector_names = v_names;
        tensor_names = t_names;
        surface_names = surf_names;

        if (!xcg)     throw std::invalid_argument("Visualization: xcg pointer is null.");
        if (!cgelcon) throw std::invalid_argument("Visualization: cgelcon pointer is null.");
        if (!telem)   throw std::invalid_argument("Visualization: telem pointer is null.");
        if (nd != 2 && nd != 3) throw std::invalid_argument("nd must be 2 or 3.");
        if (npoints < 0 || ne <= 0 || npe <= 0 || nce <= 0 || nve <= 0)
            throw std::invalid_argument("Visualization: invalid sizes (np, ne, npe, nce, nve).");

        // cgnodes padded to 3D
        cgnodes.assign(3 * npoints, 0.0f);
        for (int p = 0; p < npoints; ++p) {
            const int in_col  = p * nd;
            const int out_col = p * 3;
            cgnodes[out_col + 0] = (float) xcg[in_col + 0];
            cgnodes[out_col + 1] = (nd >= 2) ? (float) xcg[in_col + 1] : 0.0f;
            cgnodes[out_col + 2] = (nd == 3) ? (float) xcg[in_col + 2] : 0.0f;
        }

        // connectivity
        ncells = nce * ne;
        cgcells.assign(nve * ncells, 0);
        const int plane = nve * nce;
        for (int el = 0; el < ne; ++el)
            for (int j2 = 0; j2 < nce; ++j2)
                for (int i2 = 0; i2 < nve; ++i2) {
                    const int r   = telem[j2 + i2 * nce];
                    const int val = cgelcon[r + el * npe];
                    cgcells[i2 + j2 * nve + el * plane] = static_cast<int32_t>(val);
                }

        // offsets & types
        celloffsets.resize(ncells);
        for (int c = 0; c < ncells; ++c) celloffsets[c] = (c + 1) * nve;

        const int celltype = (nd == 2) ? ((elemtype == 0) ? 5 : 9)  // tri/quad
                                       : ((elemtype == 0) ? 10 : 12); // tet/hex
        celltypes.assign(ncells, static_cast<uint8_t>(celltype));

        // appended offsets
        const int fbytesize = 4; // Float32
        const int ibytesize = 4; // Int32
        const int obytesize = 8; // UInt64

        std::uint64_t offset = 0;
        auto add_off = [&](int payload_bytes) {
            std::uint64_t here = offset;
            offset += (std::uint64_t)payload_bytes + (std::uint64_t)obytesize;
            return here;
        };

        scalar_offsets.clear();
        vector_offsets.clear();
        tensor_offsets.clear();
        scalar_offsets.reserve((int)scalar_names.size());
        vector_offsets.reserve((int)vector_names.size());
        tensor_offsets.reserve((int)tensor_names.size());

        for (int i = 0; i < (int)scalar_names.size(); ++i)
            scalar_offsets.push_back(add_off(npoints * fbytesize));
        for (int i = 0; i < (int)vector_names.size(); ++i)
            vector_offsets.push_back(add_off(3 * npoints * fbytesize));
        for (int i = 0; i < (int)tensor_names.size(); ++i)
            tensor_offsets.push_back(add_off(ntc * npoints * fbytesize));

        points_offset = add_off(3 * npoints * fbytesize);
        conn_offset   = add_off(ncells * nve * ibytesize);
        offs_offset   = add_off(ncells * ibytesize);
        types_offset  = add_off(ncells * 1);
    }

    // endianness
    static const char* vtk_byte_order() {
        const uint16_t x = 0x0102;
        return (*reinterpret_cast<const uint8_t*>(&x) == 0x02)
            ? "LittleEndian" : "BigEndian";
    }

    // VTK appended block (UInt64 length + payload)
    static void write_block(std::ofstream& s, const void* data, int nbytes) {
        std::uint64_t nb = static_cast<std::uint64_t>(nbytes);
        s.write(reinterpret_cast<const char*>(&nb), sizeof(std::uint64_t));
        if (nbytes) s.write(reinterpret_cast<const char*>(data), nbytes);
        if (!s) throw std::runtime_error("Error writing appended block.");
    }

    // Write N values as Float32 appended block from either float* or double*
    template <class T>
    static void write_as_float32(std::ofstream& os, const T* src, int N) {
        using U = std::remove_cv_t<T>;
        if constexpr (std::is_same_v<U, float>) {
            // Fast path: already Float32, write directly
            write_block(os, src, N * (int)sizeof(float));
        } else {
            // U is double: convert once into a reusable buffer and write
            static thread_local std::vector<float> buf; // reuse to avoid re-allocs
            buf.resize(N);
            std::transform(src, src + N, buf.begin(),
                           [](double v){ return static_cast<float>(v); });
            write_block(os, buf.data(), N * (int)sizeof(float));
        }
    }
    
    // utilities
    static std::string rank_tag(int rank, int width = 5) {
        std::ostringstream ss; ss << "_" << std::setw(width) << std::setfill('0') << rank; return ss.str();
    }
    static std::string step_tag(int k, int width = 6) {
        std::ostringstream ss; ss << std::setw(width) << std::setfill('0') << k; return ss.str();
    }

    // PVTU: advertise ntc components for tensors
    static void write_pvtu(const std::string& pvtu_basename_no_ext,
                           const std::vector<std::string>& piece_files,
                           const std::vector<std::string>& scalar_names,
                           const std::vector<std::string>& vector_names,
                           const std::vector<std::string>& tensor_names,
                           int ntc /* nd*nd */)
    {
        const std::string fname = pvtu_basename_no_ext + ".pvtu";
        std::ofstream os(fname, std::ios::binary);
        if (!os) throw std::runtime_error("Cannot open PVTU file: " + fname);

        os << "<?xml version=\"1.0\"?>\n";
        os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\""
           << vtk_byte_order() << "\" header_type=\"UInt64\">\n";
        os << "  <PUnstructuredGrid GhostLevel=\"0\">\n";

        if (!scalar_names.empty() || !vector_names.empty() || !tensor_names.empty()) {
            os << "    <PPointData Scalars=\"scalars\">\n";
            for (const auto& s : scalar_names)
                os << "      <PDataArray type=\"Float32\" Name=\"" << s
                   << "\" NumberOfComponents=\"1\"/>\n";
            for (const auto& v : vector_names)
                os << "      <PDataArray type=\"Float32\" Name=\"" << v
                   << "\" NumberOfComponents=\"3\"/>\n";
            for (const auto& t : tensor_names)
                os << "      <PDataArray type=\"Float32\" Name=\"" << t
                   << "\" NumberOfComponents=\"" << ntc << "\"/>\n";
            os << "    </PPointData>\n";
        }

        os << "    <PPoints>\n";
        os << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
        os << "    </PPoints>\n";

        os << "    <PCells>\n";
        os << "      <PDataArray type=\"Int32\"  Name=\"connectivity\"/>\n";
        os << "      <PDataArray type=\"Int32\"  Name=\"offsets\"/>\n";
        os << "      <PDataArray type=\"UInt8\"  Name=\"types\"/>\n";
        os << "    </PCells>\n";

        for (const auto& pf : piece_files)
            os << "    <Piece Source=\"" << pf << "\"/>\n";

        os << "  </PUnstructuredGrid>\n";
        os << "</VTKFile>\n";
        if (!os) throw std::runtime_error("Error writing PVTU file.");
    }
};

