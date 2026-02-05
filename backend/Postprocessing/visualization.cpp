#include <numeric>

#ifdef USE_FLOAT
typedef float dstype;
#else
typedef double dstype; //  double is default precision 
#endif

using namespace std;


template <typename T> static void TemplateMalloc(T **data, int n)
{    
    if (n>0) *data = (T *) malloc(n*sizeof(T));       
}

template <typename T> static void TemplateFree(T *data)
{
    if (data != nullptr) {                                                      
        free(data);                                                          
        data = nullptr;                                                         
    }                                                                     
}

std::string make_path(const std::string& str1, const std::string& str2) {
    std::filesystem::path base = str1;
    std::filesystem::path tail = str2;

    // If tail is absolute, strip its root so it becomes relative    
    tail = tail.relative_path();

    std::filesystem::path full = base / tail;
    return full.lexically_normal().string();
}

int index4D(int i, int j, int k, int l, const vector<int>& shape) {
    // Column-major indexing: idx = i + j*n1 + k*n1*n2 + l*n1*n2*n3
    return i + shape[0] * (j + shape[1] * (k + shape[2] * l));
}

void masternodes(vector<dstype>& pelem, vector<int>& telem,
                 vector<dstype>& pface, vector<int>& tface,
                 vector<int>& perm, int porder, int dim, int elemtype, const std::string filename) 
{
    
    ifstream file(filename, ios::binary);
    
    if (!file) error("Error opening file: " + filename);

    // Read the full file into a vector
    file.seekg(0, ios::end);
    size_t num_bytes = file.tellg();
    file.seekg(0, ios::beg);
    size_t num_doubles = num_bytes / sizeof(dstype);

    vector<dstype> tmp(num_doubles);
    file.read(reinterpret_cast<char*>(tmp.data()), num_bytes);
    file.close();

    // Parse header
    int ndims = static_cast<int>(tmp[0]);  
    
    vector<int> narrays(ndims);
    for (int i = 0; i < ndims; ++i)
        narrays[i] = static_cast<int>(tmp[1 + i]);
    
    int offset = 1 + ndims;
    int total_blocks = 1;
    for (int d : narrays)
        total_blocks *= d;

    vector<int> sz1(total_blocks), sz2(total_blocks);
    for (int i = 0; i < total_blocks; ++i)
        sz1[i] = static_cast<int>(tmp[offset + i]);
    for (int i = 0; i < total_blocks; ++i)
        sz2[i] = static_cast<int>(tmp[offset + total_blocks + i]);

    vector<int> sz(total_blocks);
    for (int i = 0; i < total_blocks; ++i)
        sz[i] = sz1[i] * sz2[i];
    
    // cumulative offsets
    vector<int> lz(total_blocks + 1, 0);
    partial_sum(sz.begin(), sz.end(), lz.begin() + 1);
    
    // Starting point of real data
    int data_start = offset + 2 * total_blocks;
    
    auto extract_block = [&](int i, vector<dstype>& out) {
        int e = elemtype + 1;        
        int idx = index4D(i, e - 1, porder - 1, dim - 1, narrays);
        int start = lz[idx];
        int count = sz[idx];
        //printf("i = %d, e = %d, porder = %d, dim = %d, idx = %d, start = %d, count = %d\n", i, e, porder, dim, idx, start, count);
        out.resize(count);
        copy(tmp.begin() + data_start + start,
             tmp.begin() + data_start + start + count,
             out.begin());
    };

    vector<dstype>  telemd, tfaced, permd;
    
    extract_block(0, pelem);
    extract_block(1, telemd);
    extract_block(2, pface);
    extract_block(3, tfaced);
    extract_block(4, permd);
    
    if (dim==1) {
      pface.resize(1); pface[0] = 0;
      tfaced.resize(1); tfaced[0] = 1;
    }
    
    perm.resize(permd.size());    
    telem.resize(telemd.size());    
    tface.resize(tfaced.size());    
    
    for (int i=0; i<permd.size(); i++) perm[i] = (int) permd[i]-1;     
    for (int i=0; i<telemd.size(); i++) telem[i] = (int) telemd[i]-1;           
    for (int i=0; i<tfaced.size(); i++) tface[i] = (int) tfaced[i]-1;                
}

static inline int equal_row_tol(const double* a, const double* b, int dim, double tol) {
    for (int i = 0; i < dim; ++i)
        if (std::fabs(a[i] - b[i]) > tol) return 0;
    return 1;
}

static inline void fetch_node(double* out,
                              const double* dgnodes,
                              int npe, int dim, int ne,  // dims
                              int e, int a)              // element, local node
{
    // dgnodes is [npe][dim][ne] in column-major as in your routine
    for (int d = 0; d < dim; ++d)
        out[d] = dgnodes[a + npe * (d + dim * e)];
}

// Key for a grid cell: integer coordinates per dimension
struct CellKey {
    std::array<int64_t,3> c{{0,0,0}};
    uint8_t dim{0};
    bool operator==(const CellKey& o) const noexcept {
        if (dim != o.dim) return false;
        for (uint8_t i = 0; i < dim; ++i) if (c[i] != o.c[i]) return false;
        return true;
    }
};
struct CellKeyHash {
    size_t operator()(const CellKey& k) const noexcept {
        // 64-bit mix of 1..3 int64 coordinates + dim
        uint64_t h = 0x9e3779b97f4a7c15ull ^ k.dim;
        auto mix = [](uint64_t x)->uint64_t{
            x += 0x9e3779b97f4a7c15ull;
            x = (x ^ (x>>30)) * 0xbf58476d1ce4e5b9ull;
            x = (x ^ (x>>27)) * 0x94d049bb133111ebull;
            return x ^ (x>>31);
        };
        for (uint8_t i = 0; i < k.dim; ++i) {
            h ^= mix(static_cast<uint64_t>(k.c[i])) + (h<<6) + (h>>2);
        }
        return static_cast<size_t>(h);
    }
};

// Build a cell key from a point by snapping to a grid of size h (≈ tol)
static inline CellKey make_cell_key(const double* p, int dim, double inv_h) {
    CellKey key; key.dim = static_cast<uint8_t>(dim);
    // Use floor(p/h) to avoid rounding-edge inconsistencies; we will still check neighbors.
    for (int d = 0; d < dim; ++d) {
        // avoid overflow by using floor(p * inv_h)
        double s = std::floor(p[d] * inv_h);
        key.c[d] = static_cast<int64_t>(s);
    }
    return key;
}

int mkelconcg_hashgrid(
    double* cgnodes,       // [dim][max_nodes], user-allocated
    int*    cgelcon,       // [npe][ne], user-allocated
    const double* dgnodes, // [npe][dim][ne]
    int npe, int dim, int ne,
    double tol = 1e-8)
{
    // Grid cell size. Using h = tol ensures any two points within tol
    // are in the same or neighboring cells (so we check neighbors).
    const double h = tol;
    const double inv_h = 1.0 / h;

    // cell -> list of CG node indices (candidates for matching)
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> grid;
    // Reserve roughly number of points to reduce rehashing
    grid.reserve(static_cast<size_t>(ne) * static_cast<size_t>(npe));

    int ncg = 0;
    double p[3]; // up to 3D

    // Neighbor cell offsets: all combinations in {-1,0,1}^dim
    std::array<int,27> off; // max 3^3
    // We'll iterate via nested loops per dim, but generate on the fly.

    for (int e = 0; e < ne; ++e) {
        for (int a = 0; a < npe; ++a) {
            fetch_node(p, dgnodes, npe, dim, ne, e, a);

            // current cell
            CellKey ck = make_cell_key(p, dim, inv_h);

            int found = -1;

            // Visit neighbor cells (3^dim):
            if (dim == 1) {
                for (int dx = -1; dx <= 1 && found < 0; ++dx) {
                    CellKey q = ck; q.c[0] += dx;
                    auto it = grid.find(q);
                    if (it == grid.end()) continue;
                    for (int idx : it->second) {
                        if (equal_row_tol(&cgnodes[idx * dim], p, dim, tol)) { found = idx; break; }
                    }
                }
            } else if (dim == 2) {
                for (int dx = -1; dx <= 1 && found < 0; ++dx) {
                    for (int dy = -1; dy <= 1 && found < 0; ++dy) {
                        CellKey q = ck; q.c[0] += dx; q.c[1] += dy;
                        auto it = grid.find(q);
                        if (it == grid.end()) continue;
                        for (int idx : it->second) {
                            if (equal_row_tol(&cgnodes[idx * dim], p, dim, tol)) { found = idx; break; }
                        }
                    }
                }
            } else { // dim == 3
                for (int dx = -1; dx <= 1 && found < 0; ++dx) {
                    for (int dy = -1; dy <= 1 && found < 0; ++dy) {
                        for (int dz = -1; dz <= 1 && found < 0; ++dz) {
                            CellKey q = ck; q.c[0] += dx; q.c[1] += dy; q.c[2] += dz;
                            auto it = grid.find(q);
                            if (it == grid.end()) continue;
                            for (int idx : it->second) {
                                if (equal_row_tol(&cgnodes[idx * dim], p, dim, tol)) { found = idx; break; }
                            }
                        }
                    }
                }
            }

            const int aidx = a + npe * e;
            if (found >= 0) {
                // reuse existing CG node
                cgelcon[aidx] = found;
            } else {
                // create new CG node
                for (int d = 0; d < dim; ++d) cgnodes[ncg * dim + d] = p[d];
                cgelcon[aidx] = ncg;

                // insert into its cell’s list
                auto& vec = grid[ck];
                vec.push_back(ncg);

                ++ncg;
            }
        }
    }

    return ncg;
}

int mkent2elem(
    int*& rowent2elem, // output [ndof+1], allocated inside
    int*& colent2elem,  // output [nnz], allocated inside        
    const int* cgelcon, // [nrow * ne], column-major
    int nrow,
    int ne
) {
    int total = nrow * ne;
    int entmax = 0;
    for (int i = 0; i < total; ++i)
        if (cgelcon[i] > entmax)
            entmax = cgelcon[i];

    bool* mark = (bool*)calloc(entmax + 1, sizeof(bool));
    for (int i = 0; i < total; ++i)
        if (cgelcon[i] >= 0)
            mark[cgelcon[i]] = true;

    int* ent2ind = (int*)malloc((entmax + 1) * sizeof(int));
    int ndof = 0;
    for (int i = 0; i <= entmax; ++i) {
        if (mark[i]) {
            ent2ind[i] = ndof++;
        } else {
            ent2ind[i] = -1;
        }
    }

    TemplateMalloc(&rowent2elem, (ndof+1));
    for (int i=0; i<ndof+1; i++) rowent2elem[i]=0;

    int* counter = (int*)calloc(ndof, sizeof(int));
    bool* seen = (bool*)calloc(entmax + 1, sizeof(bool));

    for (int e = 0; e < ne; ++e) {
        for (int i = 0; i < nrow; ++i) {
            int ent = cgelcon[i + nrow * e];
            if (ent >= 0 && !seen[ent]) {
                seen[ent] = true;
                rowent2elem[ent2ind[ent] + 1]++;
            }
        }
        for (int i = 0; i < nrow; ++i) {
            int ent = cgelcon[i + nrow * e];
            if (ent >= 0)
                seen[ent] = false;
        }
    }

    for (int i = 1; i <= ndof; ++i)
        rowent2elem[i] += rowent2elem[i - 1];

    TemplateMalloc(&colent2elem, rowent2elem[ndof]);
    for (int i=0; i<rowent2elem[ndof]; i++) colent2elem[i]=0;

    for (int e = 0; e < ne; ++e) {
        for (int i = 0; i < nrow; ++i) {
            int ent = cgelcon[i + nrow * e];
            if (ent >= 0 && !seen[ent]) {
                seen[ent] = true;
                int id = ent2ind[ent];
                colent2elem[rowent2elem[id] + counter[id]++] = e;
            }
        }
        for (int i = 0; i < nrow; ++i) {
            int ent = cgelcon[i + nrow * e];
            if (ent >= 0)
                seen[ent] = false;
        }
    }

    TemplateFree(mark);
    TemplateFree(counter);
    TemplateFree(seen);
    TemplateFree(ent2ind);
    
    return ndof;
}

int xiny(const double* xcg, const double* xdg, int npe, int dim) 
{
    for (int i = 0; i < npe; ++i) {
        int match = 1;
        for (int d = 0; d < dim; ++d) {
            if (fabs(xcg[d] - xdg[i + npe * d]) > 1e-8) {
                match = 0;
                break;
            }
        }
        if (match) return i;
    }
    return -1; // not found
}

void map_cgent2dgent(
    int*& cgent2dgent,        // [nnz], output mapping (same shape as colent2elem)
    const int* rowent2elem, // [nent+1]
    const int* colent2elem, // [nnz]    
    const double* cgnodes,  // [nent * dim], row-major
    const double* dgnodes,  // [npe * dim * ne], column-major
    int npe, int dim, int nent) 
{    
    TemplateMalloc(&cgent2dgent, rowent2elem[nent]);
    for (int i = 0; i < nent; ++i) {
        int start = rowent2elem[i];
        int end = rowent2elem[i + 1];
        const double* xcg = &cgnodes[i * dim];
        for (int j = start; j < end; ++j) {
            int elem = colent2elem[j];
            const double* xdg = &dgnodes[npe * dim * elem];
            int in = xiny(xcg, xdg, npe, dim);
            if (in==-1) {
              for (int k = 0; k < npe; ++k) 
                for (int d = 0; d < dim; ++d) 
                  printf("%d %d: %g\n", k, d, fabs(xcg[d] - xdg[k + npe * d]));
              error("xcg and xdg do not match.");            
            }        
            cgent2dgent[j] = elem * npe + in; // global DG index
        }
    }
}

void GetElemNodes(dstype* unView, const dstype* uView, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
{
    int nn = np * (e2 - e1);
    int ncu = nc2 - nc1;
    int N = nn * ncu;
    int K = np * nc;
    
    for (int idx=0; idx < N; idx++) {
        int i = idx % nn;  // [0, np*(e2-e1)]
        int j = idx / nn;  // [0, ncu]
        int k = i % np;    // [0, np]
        int e = i / np + e1;
        unView[idx] = uView[k + (j + nc1) * np + e * K];
    }
}

void VisDG2CG(float* ucg, const dstype* udg, const int* cgent2dgent, const int* colent2elem, const int* rowent2elem, int ne1, const int ncg, int ndg, int na, int nb, int nc)
{        
    int N = nb * ncg * nc;    
    for (int idx=0; idx < N; idx++) {
        int i1 = idx%nb; 
        int n  = idx/nb; 
        int i  = n%ncg;
        int i3 = n/ncg;
        int ic = i1 + nb*i3;
        // idx = i1 + nb*i + nb*ncg*i3

        dstype sum = 0.0, n1 = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];               
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            if (e<ne1) {
                sum += udg[cgent2dgent[rowent2elem[i]+k] + ndg*ic]; 
                n1 += 1.0;
            }
        }
        ucg[i1 + na*i + na*ncg*i3] = (float) (sum/n1);
    }
}

    // common.nc = app.ndims[5]; // number of compoments of (u, q)
    // common.ncu = app.ndims[6];// number of compoments of (u)        
    // common.ncq = app.ndims[7];// number of compoments of (q)
    // //common.ncp = app.ndims[8];// number of compoments of (p)    
    // common.nco = app.ndims[9];// number of compoments of (o)    
    // common.nch = app.ndims[10];// number of compoments of (uhat)
    // common.ncx = app.ndims[11];// number of compoments of (xdg)        
    // common.nce = app.ndims[12];// number of compoments of (output)        
    // common.ncw = app.ndims[13];//number of compoments of (w)
    // common.nsca = app.ndims[14];// number of components of scalar fields for visualization
    // common.nvec = app.ndims[15];// number of components of vector fields for visualization
    // common.nten = app.ndims[16];// number of components of tensor fields for visualization
    // common.nsurf = app.ndims[17];// number of components of surface fields for visualization, storage, and QoIs
    // common.nvqoi = app.ndims[18];// number of volume quantities of interest (QoIs)    

    // common.nd = master.ndims[0];     // spatial dimension    
    // common.elemtype = master.ndims[1]; 
    // common.nodetype = master.ndims[2]; 
    // common.porder = master.ndims[3]; 
    // common.pgauss = master.ndims[4]; 
    // common.npe = master.ndims[5]; // number of nodes on master element
    // common.npf = master.ndims[6]; // number of nodes on master face       
    // common.nge = master.ndims[7]; // number of gauss points on master element
    // common.ngf = master.ndims[8]; // number of gauss poInts on master face          
    // common.np1d = master.ndims[9]; // number of node points on 1D element
    // common.ng1d = master.ndims[10]; // number of gauss poInts on 1D element          

struct appsolstruct {     
    string exasimpath = "";  
    string filein = "";       // Name of binary file with input data
    string fileout = "";      // Name of binary file to write the solution                
    int modelnumber;
    int mpiRank = 0;
    int mpiProcs = 0;
    int saveSolFreq = 0;
    int tdep = 0;
    int currentstep = 0;
    int tsteps = 0;
    int timestepOffset = 0;
    int nd = 0;
    int nc = 0;
    int ncu = 0;
    int ncw = 0;
    int ncx = 0;
    int nco = 0;
    int szxcg = 0;
    int porder = 0;
    int nsca = 0;
    int nvec = 0;
    int nten = 0;
    int nsurf = 0;
    int npe = 0;
    int ne1 = 0;
    int elemtype = 0;          
    dstype time;
    dstype* dt;
    dstype* physicsparam;
    dstype* externalparam;
    dstype *xdg=nullptr; // spatial coordinates
    dstype *udg=nullptr; // solution (u, q) 
    dstype *odg=nullptr; // auxilary term 
    dstype *wdg=nullptr; // dw/dt = u (wave problem)  
};

class CVisualization {
public:
    dstype* xdg=nullptr;
    dstype* udg=nullptr;
    dstype* vdg=nullptr;
    dstype* wdg=nullptr;
    dstype* fdg=nullptr;

    float* scafields=nullptr;
    float* vecfields=nullptr;
    float* tenfields=nullptr;
    float* srffields=nullptr;

    int *rowent2elem=nullptr;
    int *cgent2dgent=nullptr;
    int *colent2elem=nullptr;

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
    CVisualization(appsolstruct& appsol) {      
        rank = appsol.mpiRank;
        int nd_in   = appsol.nd;
        int npoints_in = appsol.szxcg / nd_in;
        
        if (npoints_in > 0 && nd_in > 1) {            
            int porder  = appsol.porder;        
            int nsca    = appsol.nsca;
            int nvec    = appsol.nvec;            
            int nten    = appsol.nten;            
            int nsurf   = appsol.nsurf;            
            int npe     = appsol.npe;
            int ne      = appsol.ne1;
            int elemtype= appsol.elemtype;
            int nve_in  = (elemtype==0) ? (nd_in + 1) : std::pow(2, nd_in);
                
            std::string fn1 = make_path(appsol.exasimpath, "text2code/text2code/masternodes.bin");
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
                                                

            int dim = nd_in;
            int ncf = max(max(nsca, 3*nvec), dim*dim*nten);          
            TemplateMalloc(&xdg, npe*ne*dim);
            TemplateMalloc(&udg, npe*ne*appsol.nc);
            TemplateMalloc(&vdg, npe*ne*appsol.nco);
            TemplateMalloc(&wdg, npe*ne*appsol.ncw);            
            TemplateMalloc(&fdg, npe*ne*ncf);

            GetElemNodes(xdg, appsol.xdg, npe, dim, 0, dim, 0, ne);
          
            dstype *xcg=nullptr;
            int *cgelcon=nullptr;
            TemplateMalloc(&xcg, npe*dim*ne);
            TemplateMalloc(&cgelcon, npe*ne);
            int ncgnodes = mkelconcg_hashgrid(xcg, cgelcon, appsol.xdg, npe, dim, ne);
            int ncgdof = mkent2elem(rowent2elem, colent2elem, cgelcon, npe, ne);
            map_cgent2dgent(cgent2dgent, rowent2elem, colent2elem, xcg, appsol.xdg, npe, dim, ncgdof);                                          
          
            Init(xcg, nd_in, npoints_in, cgelcon, npe, ne,
                 telem.data(), nce, nve_in, elemtype,
                 scalars, vectors, tensors, surfaces);

            savemode = (nsca + nvec + nten > 0); 
        
            scafields = (float *) malloc(npoints*nsca*sizeof(float));
            vecfields = (float *) malloc(3*npoints*nvec*sizeof(float));
            tenfields = (float *) malloc(ntc*npoints*nten*sizeof(float));
            host_alloc_backend = 0;
            
            TemplateFree(xcg);
            TemplateFree(cgelcon);

            //cout<<ne<<"  "<<npoints<<endl;
            if (appsol.mpiRank == 0) printf("finish CVisualization constructor... \n");    
        }        
    }

    ~CVisualization() {
        // Free scafields / vecfields / tenfields / srffields according to how they were allocated.
        TemplateFree(scafields);
        TemplateFree(vecfields);
        TemplateFree(tenfields);
        TemplateFree(srffields); 
        TemplateFree(xdg);
        TemplateFree(udg);
        TemplateFree(vdg);
        TemplateFree(wdg);
        TemplateFree(wdg);
        TemplateFree(cgent2dgent);
        TemplateFree(rowent2elem);
        TemplateFree(colent2elem);
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

        if (!xcg)     error("Visualization: xcg pointer is null.");
        if (!cgelcon) error("Visualization: cgelcon pointer is null.");
        if (!telem)   error("Visualization: telem pointer is null.");
        if (nd != 2 && nd != 3) throw std::invalid_argument("nd must be 2 or 3.");
        if (npoints < 0 || ne <= 0 || npe <= 0 || nce <= 0 || nve <= 0)
            error("Visualization: invalid sizes (np, ne, npe, nce, nve).");

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

void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
}

void KokkosVisTensors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
}

void KokkosVisVectors(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
}

void WriteParaview(CVisualization& vis, appsolstruct& appsol) 
{
    // Decide whether we should write a file on this step
    bool writeSolution = false;
    
    if (appsol.tdep == 1) {
       if (appsol.currentstep==0 && appsol.mpiRank==0) {
          string ext = (appsol.mpiProcs==1) ? "vtu" : "pvtu";                                  
          vis.pvdwrite_series(appsol.fileout + "vis", appsol.dt, appsol.tsteps, appsol.saveSolFreq, ext);                          
       }
        
        // Time-dependent: only write every 'saveSolFreq' steps
        writeSolution = ((appsol.currentstep + 1) % appsol.saveSolFreq) == 0;        
    } else {
        // Steady / not time-dependent: always write
        writeSolution = true;
    }

   if (writeSolution) { 
       int nd = appsol.nd;   
       int nc = appsol.nc;  
       int ncu = appsol.ncu;  
       int ncx = appsol.ncx;        
       int nco = appsol.nco;  
       int ncw = appsol.ncw;  
       int nsca = appsol.nsca; 
       int nvec = appsol.nvec;  
       int nten = appsol.nten;     
       int npe  = appsol.npe;     
       int ne   = appsol.ne1;      
       int ndg  = npe * ne;
       int ncg  = vis.npoints;
           
       GetElemNodes(vis.udg, appsol.udg, npe, nc, 0, nc, 0, ne);
       if (nco > 0) GetElemNodes(vis.vdg, appsol.odg, npe, nco, 0, nco, 0, ne);
       if (ncw > 0) GetElemNodes(vis.wdg, appsol.wdg, npe, ncw, 0, ncw, 0, ne);
       
       int numPoints = npe*ne;
       if (nsca > 0) {        
            KokkosVisScalars(vis.fdg, vis.xdg, vis.udg, vis.vdg, vis.wdg, appsol.externalparam, appsol.physicsparam, 
                       appsol.time, appsol.modelnumber, npe*ne, nc, ncu, nd, ncx, nco, ncw);                       
            VisDG2CG(vis.scafields, vis.fdg, vis.cgent2dgent, vis.colent2elem, vis.rowent2elem, ne, ncg, ndg, 1, 1, nsca);
       }    
       if (nvec > 0) {        
            KokkosVisVectors(vis.fdg, vis.xdg, vis.udg, vis.vdg, vis.wdg, appsol.externalparam, appsol.physicsparam, 
                       appsol.time, appsol.modelnumber, npe*ne, nc, ncu, nd, ncx, nco, ncw);                       
            VisDG2CG(vis.vecfields, vis.fdg, vis.cgent2dgent, vis.colent2elem, vis.rowent2elem, ne, ncg, ndg, 3, ncx, nvec);
       }
       if (nten > 0) {        
            KokkosVisTensors(vis.fdg, vis.xdg, vis.udg, vis.vdg, vis.wdg, appsol.externalparam, appsol.physicsparam, 
                       appsol.time, appsol.modelnumber, npe*ne, nc, ncu, nd, ncx, nco, ncw);                       
            VisDG2CG(vis.tenfields, vis.fdg, vis.cgent2dgent, vis.colent2elem, vis.rowent2elem, ne, ncg, ndg, vis.ntc, vis.ntc, nvec);
       }

       string baseName = appsol.fileout + "vis";
       if (appsol.tdep == 1) {
           std::ostringstream ss; 
           ss << std::setw(6) << std::setfill('0') << appsol.currentstep+appsol.timestepOffset+1; 
           baseName = baseName + "_" + ss.str();           
       }
       
       if (appsol.mpiProcs==1)                
            vis.vtuwrite(baseName, vis.scafields, vis.vecfields, vis.tenfields);
       else
            vis.vtuwrite_parallel(baseName, appsol.mpiRank, appsol.mpiProcs, vis.scafields, vis.vecfields, vis.tenfields);
   }
}
