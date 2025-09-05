/**
 * @file connectivity.cpp
 * @brief Mesh connectivity routines for finite element/DG solvers.
 *
 * This file provides functions and structures for building and managing mesh connectivity,
 * including element-to-node, face-to-element, and boundary condition mappings.
 *
 * Main components:
 * - Conn struct: Holds all connectivity arrays and metadata for a mesh partition.
 * - apply_bcm: Applies boundary condition map to mesh faces.
 * - fix_f2t_consistency: Ensures face-to-element mapping consistency for parallel partitions.
 * - build_connectivity: Constructs element and face connectivity arrays for DG/CG methods.
 * - connectivity: Main routine to partition faces, apply boundary conditions, and build connectivity.
 * - mkelconcg: Generates unique CG node coordinates and element-to-node mapping.
 * - mkent2elem: Builds node-to-element connectivity in CRS format.
 * - map_cgent2dgent: Maps CG nodes to DG nodes using CRS connectivity.
 * - mkdge2dgf: Constructs CRS mapping from face connectivity to per-entity DOFs.
 * - removeBoundaryFaces: Removes boundary faces from face connectivity based on block info.
 * - divide_interval, mkfaceblocks: Utility routines for partitioning elements/faces into blocks.
 * - buildConn: High-level routine to build all mesh connectivity for a PDE, mesh, and master element.
 *
 * Dependencies:
 * - Requires PDE, Mesh, Master, DMD structures (not defined here).
 * - Relies on utility functions: select_columns, find, permute_columns, simple_bubble_sort,
 *   unique_count, unique_ints, print2iarray, print2darray, error, CPUFREE, etc.
 *
 * Usage:
 * - Call buildConn to initialize a Conn object for a given mesh and PDE setup.
 * - Use connectivity and related routines for custom mesh partitioning and connectivity construction.
 *
 * Notes:
 * - Supports both serial and parallel (MPI) mesh partitions.
 * - Handles both DG and CG node mappings.
 * - Designed for high-order finite element/DG solvers.
 */

#ifndef __CONNECTIVITY
#define __CONNECTIVITY
    
void select_columns(int* a_new, const int* a, const int* ind, int m, int k) 
{
    for (int j = 0; j < k; ++j) {
        int col = ind[j];
        for (int i = 0; i < m; ++i) {
            a_new[i + j * m] = a[i + col * m];
        }
    }
}

int find(const int* a, int b, int m, int n, int k, int opts) 
{
    int count = 0;
    if (opts==0) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] == b) count++;                        
      }
    }
    else if (opts==1) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] <= b) count++;
      }
    }
    else if (opts==2) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] >= b) count++;
      }
    }    
    return count;
}

int find(int* indices, const int* a, int b, int m, int n, int k, int opts) 
{
    int count = 0;
    if (opts==0) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] == b) {
              indices[count++] = i;
          }
      }
    }
    else if (opts==1) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] <= b) {
              indices[count++] = i;
          }
      }
    }
    else if (opts==2) {
      for (int i = 0; i < n; ++i) {
          if (a[k + i*m] >= b) {
              indices[count++] = i;
          }
      }
    }
    
    return count;
}

int unique_count(int *b, int *c, const int *a, int n)
{
    if (n == 0) return 0;

    int uniq = 0;          /* index in b/c */
    double current = a[0];
    int     count   = 1;

    for (int i = 1; i < n; ++i) {
        if (a[i] == current) {
            ++count;          /* same value – just bump counter */
        } else {
            /* flush previous run */
            b[uniq] = current;
            c[uniq] = count;
            ++uniq;

            /* start new run */
            current = a[i];
            count   = 1;
        }
    }
    /* flush final run */
    b[uniq] = current;
    c[uniq] = count;
    ++uniq;

    return uniq;              /* number of unique values */
}

void permute_columns(int* a, const int* ind, int m, int k) 
{    
  if ((k > 0) && (m > 0)) {
    int* a_new = (int*)malloc(m * k * sizeof(int));
    
    for (int j = 0; j < k; ++j) {
        int col = ind[j];
        for (int i = 0; i < m; ++i) {
            a_new[i + j * m] = a[i + col * m];
        }
    }
    
    for (int j = 0; j < k; ++j)
      for (int i = 0; i < m; ++i) 
        a[i + j * m] = a_new[i + j * m];
    
    free(a_new);
  }
}

int unique_ints(int* arr, int n) {
    if (n <= 1) return n;
    std::sort(arr, arr + n);                       // introsort (quick+heap+insertion)
    int* last = std::unique(arr, arr + n);         // compacts uniques to front
    return int(last - arr);                        // new length m
}

// Hashable key of fixed max size 8 (adjust if needed)
struct FaceKey {
    std::array<int,8> a{};
    uint8_t len{}; // nnf

    bool operator==(const FaceKey& o) const noexcept {
        if (len != o.len) return false;
        for (uint8_t i = 0; i < len; ++i) if (a[i] != o.a[i]) return false;
        return true;
    }
};
struct FaceKeyHash {
    size_t operator()(const FaceKey& k) const noexcept {
        // simple 64-bit mix across len entries
        uint64_t h = 0x9E3779B97F4A7C15ull ^ k.len;
        for (uint8_t i = 0; i < k.len; ++i) {
            uint64_t x = static_cast<uint64_t>(k.a[i]) + 0x9E3779B97F4A7C15ull;
            x ^= (x >> 30); x *= 0xBF58476D1CE4E5B9ull;
            x ^= (x >> 27); x *= 0x94D049BB133111EBull;
            x ^= (x >> 31);
            h ^= x + 0x9E3779B97F4A7C15ull + (h<<6) + (h>>2);
        }
        return static_cast<size_t>(h);
    }
};

int mkf2e_hash(int* f2e,
               const int* e2n, const int* local_faces,
               int ne, int nne, int nnf, int nfe)
{
    using Val = std::pair<int,int>; // (elem, lf)
    std::unordered_map<FaceKey, Val, FaceKeyHash> map;
    map.reserve(static_cast<size_t>(ne) * nfe * 1.3);

    int out = 0;

    auto make_key = [&](int e, int lf) {
        FaceKey k; k.len = static_cast<uint8_t>(nnf);
        // gather face nodes
        for (int i = 0; i < nnf; ++i) {
            int ln = local_faces[lf*nnf + i];
            k.a[i] = e2n[e*nne + ln];
        }
        std::sort(k.a.begin(), k.a.begin() + nnf);
        return k;
    };

    for (int e = 0; e < ne; ++e) {
        for (int lf = 0; lf < nfe; ++lf) {
            FaceKey k = make_key(e, lf);
            auto it = map.find(k);
            if (it == map.end()) {
                map.emplace(std::move(k), Val{e, lf});
            } else {
                // found the interior faces
                const auto [e0, lf0] = it->second;
                f2e[4*out + 0] = e0;  f2e[4*out + 1] = lf0;
                f2e[4*out + 2] = e;   f2e[4*out + 3] = lf;
                map.erase(it);
                ++out;
            }
        }
    }

    // Remaining entries in map are boundary faces
    for (const auto& kv : map) {
        const auto [e0, lf0] = kv.second;
        f2e[4*out + 0] = e0;  f2e[4*out + 1] = lf0;
        f2e[4*out + 2] = -1;  f2e[4*out + 3] = -1;
        ++out;
    }

    return out;
}

// ---- Build e2e: for each element 'e' and local face 'lf', write neighbor element id or -1 ----
// e2e must be sized to ne * nfe (ints).
int mke2e_hash(int* e2e,
               const int* e2n,          // [ne * nne] global node ids per element
               const int* local_faces,  // [nfe * nnf] local node indices per local face
               int ne, int nne, int nnf, int nfe)
{
    using Val = std::pair<int,int>; // (elem, lf) of the first owner of the face

    // init all neighbors to -1 (assume boundary until proven otherwise)
    for (int i = 0; i < ne * nfe; ++i) e2e[i] = -1;

    // map holds faces seen once (awaiting their twin)
    std::unordered_map<FaceKey, Val, FaceKeyHash> map;
    map.reserve(static_cast<size_t>(ne) * static_cast<size_t>(nfe) * 13 / 10); // slack to keep LF low

    int out = 0;

    auto make_key = [&](int e, int lf) {
        FaceKey k; k.len = static_cast<uint8_t>(nnf);
        // gather global node ids of this face
        for (int i = 0; i < nnf; ++i) {
            const int ln = local_faces[lf * nnf + i];   // local node id on this face
            k.a[i] = e2n[e * nne + ln];                 // map to global node id
        }
        // canonicalize (orientation-independent)
        std::sort(k.a.begin(), k.a.begin() + nnf);
        return k;
    };

    for (int e = 0; e < ne; ++e) {
        for (int lf = 0; lf < nfe; ++lf) {
            FaceKey k = make_key(e, lf);
            auto it = map.find(k);

            if (it == map.end()) {
                // first time seeing this face
                map.emplace(std::move(k), Val{e, lf});
            } else {
                // found twin: set both directions and remove from map
                const auto [e0, lf0] = it->second;

                e2e[e0 * nfe + lf0] = e;   // neighbor of the first owner is current element
                e2e[e  * nfe + lf ] = e0;  // neighbor of current element is the first owner

                map.erase(it);
                ++out;
            }
        }
    }

    return out; 
}

void apply_bcm(int* bf, const int* fi, const int* bcm, int n, int nbcm) 
{
    // initialize bf to 0
    for (int i = 0; i < n; ++i)
        bf[i] = 0;

    // apply boundary condition map
    for (int j = 0; j < nbcm; ++j) {
        for (int i = 0; i < n; ++i) {
            if (fi[i] == j) {  // MATLAB is 1-based
                bf[i] = bcm[j];
            }
        }
    }
}

void fix_f2t_consistency(int* f2t, const int* elempart, int nf) 
{
    for (int n = 0; n < nf; ++n) {
        int e1 = f2t[0 + 4 * n];
        int e2 = f2t[2 + 4 * n];

        if (e2 >= 0) {  // interior face: both elements exist
            int g1 = elempart[e1];
            int g2 = elempart[e2];

            if (g2 < g1) {
                // swap rows 1:2 with 3:4 (MATLAB: f2t(3:4,n) <-> f2t(1:2,n))
                int tmp0 = f2t[0 + 4 * n];
                int tmp1 = f2t[1 + 4 * n];

                f2t[0 + 4 * n] = f2t[2 + 4 * n];
                f2t[1 + 4 * n] = f2t[3 + 4 * n];

                f2t[2 + 4 * n] = tmp0;
                f2t[3 + 4 * n] = tmp1;
            }

            // Check again
            e1 = f2t[0 + 4 * n];
            e2 = f2t[2 + 4 * n];
            g1 = elempart[e1];
            g2 = elempart[e2];

            if (g2 <= g1) {
                fprintf(stderr, "Error: DOF consistency check failed for face %d (g2 <= g1)\n", n);
                exit(1);
            }
        }
    }
}
  
void build_connectivity(int* elemcon, int* facecon, const int* f2t, const int* face, const int* t,        
                        const int* perm, const int* permind, int dim, int elemtype, int npf, int nfe, int npe, int ne, int nf) 
{
    int* facenode1 = (int*) malloc(sizeof(int) * npf);
    int* facenode2 = (int*) malloc(sizeof(int) * npf);
    for (int i = 0; i < npf; ++i) {
        facenode1[i] = i;  // zero-based indexing
        facenode2[i] = permind[i];
    }
    
    if (dim <= 2) {
        for (int i = 0; i < nf; ++i) {
            int e1 = f2t[0 + 4 * i];
            int l1 = f2t[1 + 4 * i];
            int e2 = f2t[2 + 4 * i];
            int l2 = f2t[3 + 4 * i];

            for (int j = 0; j < npf; ++j)
                facecon[0 + 2 * j + 2 * npf * i] = e1 * npe + perm[j + npf * l1];

            for (int j = 0; j < npf; ++j)
                elemcon[j + npf * (l1 + nfe * e1)] = i * npf + facenode1[j];

            if (e2 >= 0) {
                for (int j = 0; j < npf; ++j)
                    facecon[1 + 2 * j + 2 * npf * i] = e2 * npe + perm[permind[j] + npf * l2];

                for (int j = 0; j < npf; ++j)
                    elemcon[j + npf * (l2 + nfe * e2)] = i * npf + facenode2[j];
            } else {
                for (int j = 0; j < npf; ++j)
                    facecon[1 + 2 * j + 2 * npf * i] = facecon[0 + 2 * j + 2 * npf * i];
            }
        }
    } else {
        int nvf = (elemtype == 0) ? 3 : 4;
        int nve = (elemtype == 0) ? 4 : 8;      
        
        for (int i = 0; i < nf; ++i) {
            int e1 = f2t[0 + 4 * i];
            int l1 = f2t[1 + 4 * i];
            int e2 = f2t[2 + 4 * i];
            int l2 = f2t[3 + 4 * i];

            for (int j = 0; j < npf; ++j)
                facecon[0 + 2 * j + 2 * npf * i] = e1 * npe + perm[j + npf * l1];

            for (int j = 0; j < npf; ++j)
                elemcon[j + npf * (l1 + nfe * e1)] = i * npf + facenode1[j];
            
            if (e2 >= 0) {
                int f1[4], f2[4];
                for (int k = 0; k < nvf; ++k) {
                    f1[k] = t[face[k + nvf * l1] + nve * e1];
                    f2[k] = t[face[k + nvf * l2] + nve * e2];
                }

                int k = -1;
                if (elemtype == 0) {
                    // Tetrahedron face matching
                    if (f1[0] == f2[0] && f1[1] == f2[2] && f1[2] == f2[1]) k = 0;
                    else if (f1[0] == f2[1] && f1[1] == f2[0] && f1[2] == f2[2]) k = 1;
                    else if (f1[0] == f2[2] && f1[1] == f2[1] && f1[2] == f2[0]) k = 2;
                } else {
                    // Hexahedron face matching
                    if (f1[0]==f2[0] && f1[1]==f2[3] && f1[2]==f2[2] && f1[3]==f2[1]) k = 0;
                    else if (f1[0]==f2[3] && f1[1]==f2[2] && f1[2]==f2[1] && f1[3]==f2[0]) k = 1;
                    else if (f1[0]==f2[2] && f1[1]==f2[1] && f1[2]==f2[0] && f1[3]==f2[3]) k = 2;
                    else if (f1[0]==f2[1] && f1[1]==f2[0] && f1[2]==f2[3] && f1[3]==f2[2]) k = 3;
                }
                if (k < 0) {                    
                    error("Mesh connectivity is wrong\n");
                }
                for (int j = 0; j < npf; ++j)
                    facecon[1 + 2 * j + 2 * npf * i] = e2 * npe + perm[permind[j + npf * k] + npf * l2];

                for (int j = 0; j < npf; ++j)
                    elemcon[j + npf * (l2 + nfe * e2)] = i * npf + facenode1[permind[j + npf * k]];                                
            } else {
                for (int j = 0; j < npf; ++j)
                    facecon[1 + 2 * j + 2 * npf * i] = facecon[0 + 2 * j + 2 * npf * i];
            }
        }
    }

    CPUFREE(facenode1);
    CPUFREE(facenode2);
}

int connectivity(int*& elemcon, int*& facecon, int*& f2t, int*& facepartpts, int*& facepartbnd, 
                 int* sizes, const int* bf, const int* ti, const int* elempart, 
                 const int* localfaces,  const int* perm, const int* permind,  
                 const int* bcm, int dim, int elemtype, int nve, int nvf, int nfe, 
                 int npf, int npe, int ne, int ne1, int nbcm) 
{          
    int* f2t_tmp = (int*)malloc(4 * nfe * ne * sizeof(int));
    int nf1 = mkf2e_hash(f2t_tmp, ti, localfaces, ne, nve, nvf, nfe);       
    
    int* ind = (int*)malloc(nf1 * sizeof(int));
    int nf = 0;
    for (int i=0; i<nf1; i++)
    if (f2t_tmp[0 + 4*i] < ne1)
      ind[nf++] = i;                
        
    TemplateMalloc(&f2t, 4*nf, 0);
    select_columns(f2t, f2t_tmp, ind, 4, nf);        
    CPUFREE(f2t_tmp);
        
    int nifaces = find(f2t, 0, 4, nf, 2, 2);
    int nbfaces = find(f2t, -1, 4, nf, 2, 0);      
    int* ifaces = (int*)malloc(nifaces * sizeof(int));
    int* bfaces = (int*)malloc(nbfaces * sizeof(int));
    find(ifaces, f2t, 0, 4, nf, 2, 2);  // interior faces
    find(bfaces, f2t, -1, 4, nf, 2, 0); // boundary faces
            
    int* bfa = (int*)malloc(nbfaces * sizeof(int));
    for (int i=0; i<nbfaces; i++) {
        int e = f2t[0 + 4*bfaces[i]];
        int l = f2t[1 + 4*bfaces[i]];
        bfa[i] = bf[l + nfe*e];
    }                    
              
    int* sortedbfaces = (int*)malloc(nbfaces * sizeof(int));
    simple_bubble_sort(sortedbfaces, ind, bfa, nbfaces);                  
    
    int *Tfacepartbnd=nullptr, *Tfacepartpts=nullptr;
    TemplateMalloc(&Tfacepartbnd, (nbcm+1), 0);
    TemplateMalloc(&Tfacepartpts, (nbcm+1), 0);
    Tfacepartbnd[0] = 0;
    Tfacepartpts[0] = nifaces;
    int nbc = unique_count(&Tfacepartbnd[1], &Tfacepartpts[1], sortedbfaces, nbfaces);      
    //facepartpts.erase(facepartpts.begin() + 1 + nbc, facepartpts.end());
    //facepartbnd.erase(facepartbnd.begin() + 1 + nbc, facepartbnd.end());        
    TemplateMalloc(&facepartbnd, (nbc+1), 0);
    TemplateMalloc(&facepartpts, (nbc+1), 0);
    for (int i=0; i<nbc+1; i++) {
        facepartbnd[i] = Tfacepartbnd[i];
        facepartpts[i] = Tfacepartpts[i];
    }
    CPUFREE(Tfacepartbnd); CPUFREE(Tfacepartpts);

    int* bind = (int*)malloc(nf * sizeof(int));
    for (int i=0; i<nifaces; i++) bind[i] = ifaces[i];
    for (int i=0; i<nbfaces; i++) bind[nifaces+i] = bfaces[ind[i]];
    permute_columns(f2t, bind, 4, nf); 
    
    CPUFREE(ifaces); CPUFREE(bfaces); CPUFREE(bfa); CPUFREE(ind); CPUFREE(bind);  CPUFREE(sortedbfaces);                
        
    fix_f2t_consistency(f2t, elempart, nf); 
        
    TemplateMalloc(&elemcon, npf * nfe * ne, 0);
    TemplateMalloc(&facecon, npf * 2 * nf, 0);
    for (int i=0; i<npf * nfe * ne; i++)  elemcon[i] = -1;
    for (int i=0; i<npf * 2 * nf; i++)  facecon[i] = -1;
    build_connectivity(elemcon, facecon, f2t, localfaces, ti, perm, 
          permind, dim, elemtype, npf, nfe, npe, ne, nf); 
        
    sizes[0] = nf; sizes[1] = nifaces; sizes[2] = nbfaces; sizes[3] = nbc;

    return nf;      
}

int equal_row(const double* a, const double* b, int dim, double tol) {
    for (int i = 0; i < dim; ++i)
        if (fabs(a[i] - b[i]) > tol)
            return 0;
    return 1;
}

void get_node_coord(double* out, const double* dgnodes, int npe, int dim, int e, int a) {
    for (int d = 0; d < dim; ++d)
        out[d] = dgnodes[a + npe * (d + dim * e)];
}

/**
 * @param dgnodes   - [npe][dim][ne] in column-major layout
 * @param npe, dim, ne - mesh dimensions
 * @param cgnodes   - output array [dim][max_nodes], preallocated by user
 * @param cgelcon   - output array [npe][ne], preallocated by user
 * @return number of unique CG nodes (i.e., used portion of cgnodes)
 */
int mkelconcg(
    double* cgnodes,        // [dim][max_nodes], user-allocated
    int* cgelcon,            // [npe][ne], user-allocated        
    const double* dgnodes,  // [npe][dim][ne]
    int npe, int dim, int ne) 
{
    const double tol = 1e-8;
    int ncg = 0;  // number of unique CG nodes
    double node[3];  // supports up to 3D

    for (int e = 0; e < ne; ++e) {
        for (int a = 0; a < npe; ++a) {
            get_node_coord(node, dgnodes, npe, dim, e, a);

            int found = -1;
            for (int j = 0; j < ncg; ++j) {
                if (equal_row(&cgnodes[j * dim], node, dim, tol)) {
                    found = j;
                    break;
                }
            }

            int idx = a + npe * e;
            if (found >= 0) {
                cgelcon[idx] = found;
            } else {
                for (int d = 0; d < dim; ++d) cgnodes[ncg * dim + d] = node[d];
                cgelcon[idx] = ncg;
                ncg++;
            }
        }
    }

    return ncg;  // return number of unique CG nodes
}

// reuses your equal_row
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

/**
 * Spatial-hash version:
 *  - Expected O(N) (N = ne*npe) since each point checks only a tiny local neighborhood.
 *  - Exact tolerance test preserved.
 *
 * cgnodes   - [dim][max_nodes] (row-major per your use: consecutive dims for each node)
 * cgelcon   - [npe][ne] indices into cgnodes
 * returns number of unique CG nodes
 */
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

    //int* rowent2elem = (int*)calloc(ndof + 1, sizeof(int));
    //rowent2elem.resize(ndof+1);
    //std::fill(rowent2elem.begin(), rowent2elem.end(), 0);    
    TemplateMalloc(&rowent2elem, (ndof+1), 0);
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

    //colent2elem.resize(rowent2elem[ndof]);
    //std::fill(colent2elem.begin(), colent2elem.end(), 0);        
    TemplateMalloc(&colent2elem, rowent2elem[ndof], 0);
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

    CPUFREE(mark);
    CPUFREE(counter);
    CPUFREE(seen);
    CPUFREE(ent2ind);
    
    return ndof;
}

// Helper function to match a point xcg to a row in xdg: returns 0-based index
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

// Maps CG node to corresponding DG node via CRS connectivity
void map_cgent2dgent(
    int*& cgent2dgent,        // [nnz], output mapping (same shape as colent2elem)
    int*& rowent2elem, // [nent+1]
    int*& colent2elem, // [nnz]    
    const double* cgnodes,  // [nent * dim], row-major
    const double* dgnodes,  // [npe * dim * ne], column-major
    int npe, int dim, int nent) 
{    
    //cgent2dgent.resize(rowent2elem[nent]);
    //reallocate(cgent2dgent, rowent2elem[nent]);
    TemplateMalloc(&cgent2dgent, rowent2elem[nent], 0);
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

int mkdge2dgf(int*& rowdge2dgf, int*& coldge2dgf, int*& ent2ind,
              const int* facecon, int ndgf, int entmax) 
{
    // Step 1: Extract positive entries from facecon
    int* tmp = (int*)malloc(ndgf * sizeof(int));
    int count = 0;
    for (int i = 0; i < ndgf; ++i)
        if (facecon[i] >= 0)
            tmp[count++] = facecon[i];

    // Step 2: Sort and get unique entities
    int nent = unique_ints(tmp, count);

    // Step 3: Allocate and build ent2ind
    //int* ent2ind = (int*)calloc(entmax, sizeof(int));
    //ent2ind.resize(entmax);    
    //std::fill(ent2ind.begin(), ent2ind.end(), -1);
    //reallocate(ent2ind, entmax);
    TemplateMalloc(&ent2ind, entmax, 0);
    for (int i = 0; i < entmax; ++i) ent2ind[i]=-1;    
    for (int i = 0; i < nent; ++i)
        ent2ind[tmp[i]] = i;

    // Step 4: Count occurrences per entity
    //int* rowdge2dgf = (int*)calloc(nent + 1, sizeof(int));
    //rowdge2dgf.resize(nent+1);
    //std::fill(rowdge2dgf.begin(), rowdge2dgf.end(), 0);
    //reallocate(rowdge2dgf, nent+1);
    TemplateMalloc(&rowdge2dgf, nent+1, 0);
    for (int i = 0; i < nent+1; ++i) rowdge2dgf[i]=0;    
    for (int i = 0; i < ndgf; ++i) {
        int k = facecon[i];
        if (k >= 0) {
            int ind = ent2ind[k];
            rowdge2dgf[ind + 1]++;
        }
    }

    // Step 5: Compute prefix sum
    for (int i = 0; i < nent; ++i)
        rowdge2dgf[i + 1] += rowdge2dgf[i];

    // Step 6: Fill col indices
    int total = rowdge2dgf[nent];
    //int* coldge2dgf = (int*)malloc(total * sizeof(int));
    //coldge2dgf.resize(total);
    //reallocate(coldge2dgf, total);
    TemplateMalloc(&coldge2dgf, total, 0);
    int* inc = (int*)calloc(nent, sizeof(int));
    for (int i = 0; i < ndgf; ++i) {
        int k = facecon[i];
        if (k >= 0) {
            int ind = ent2ind[k];
            int pos = rowdge2dgf[ind] + inc[ind];
            coldge2dgf[pos] = i;
            inc[ind]++;
        }
    }

    CPUFREE(tmp);
    CPUFREE(inc);
    
    return nent;
}

int removeBoundaryFaces(int* facecon, const int* fblks, int m, int npf, int nf) 
{
    std::vector<bool> keep(nf, true);

    // Step 1: Identify columns to remove
    for (int i = 0; i < m; ++i) {
        int start = fblks[0 + 3*i]-1;
        int end   = fblks[1 + 3*i]-1;
        int flag  = fblks[2 + 3*i];
        if (flag > 0) {
            for (int j = start; j <= end; ++j) {
                if (j >= 0 && j < nf) {
                    keep[j] = false;
                }
            }
        }
    }

    // Step 2: Count kept columns
    int new_nf = 0;
    for (bool k : keep) {
        if (k) ++new_nf;
    }

    // Step 3: Allocate and fill new_facecon
    std::vector<int> new_facecon(npf * new_nf);
    int col_index = 0;
    for (int j = 0; j < nf; ++j) {
        if (keep[j]) {
            for (int i = 0; i < npf; ++i) {
                new_facecon[i + col_index * npf] = facecon[i + j * npf];
            }
            ++col_index;
        }
    }
    
    // Step 4: Replace facecon
    //facecon = std::move(new_facecon);    
    for (int i=0; i<(npf * new_nf); i++) facecon[i] = new_facecon[i];

    return new_nf;
}

int divide_interval(int*& intervals, int n, int m) 
{
    if (n <= 0 || m <= 0) return 0;

    int num_intervals = (n + m - 1) / m; // ceil(n/m)        
    TemplateMalloc(&intervals, 3*num_intervals, 0);
            
    int base = n / num_intervals;
    int rem = n % num_intervals;
    int current = 1;

    for (int i = 0; i < num_intervals; ++i) {
        int len = base + (i < rem ? 1 : 0);
        intervals[3 * i] = current;
        intervals[3 * i + 1] = current + len - 1;
        intervals[3 * i + 2] = 0;
        current += len;
    }

    return num_intervals;
}

int mkfaceblocks(int*& nm, const int* mf, const int* bcm, int nmf_len, int ns) 
{
    if (ns <= 0) ns = 2048;  // default value

    int max_blocks = 0;
    for (int i = 0; i < nmf_len - 1; ++i)
        max_blocks += (mf[i + 1] - mf[i] + ns - 1) / ns;
        
    TemplateMalloc(&nm, 3*max_blocks, 0);
    int count = 0;

    for (int i = 0; i < nmf_len - 1; ++i) {
        int nf = mf[i + 1] - mf[i];
        int* intervals; 
        int nblocks = divide_interval(intervals, nf, ns);

        for (int j = 0; j < nblocks; ++j) {
            int start = mf[i] + intervals[3 * j];       // 0-based
            int end   = mf[i] + intervals[3 * j + 1];    // 0-based
            nm[3 * count + 0] = start;
            nm[3 * count + 1] = end;
            nm[3 * count + 2] = bcm[i];  // boundary code
            count++;
        }
    }

    return count; 
}

/**
 * Get the local node indices for each face of a reference element.
 *
 * @param dim        Spatial dimension (1, 2, or 3)
 * @param elemtype   0 = simplex (line, tri, tet), 1 = tensor (quad, hex)
 * @param nfe     [output] number of faces
 * @param nvf    [output] number of nodes per face
 * @return           int* (column-major face matrix of size nvf × nfe)
 *                   Caller must free() the returned pointer.
 */
void getelemface(int* face, int dim, int elemtype)
{
    int nvf = (dim == 3) ? (dim + elemtype) : dim;     
    int nfe = dim + (dim-1)*elemtype + 1;    
          
    if (dim == 1) {
        //face = (int*) malloc(sizeof(int) * (*nfe) * (nvf));
        face[0] = 1;  // face 1 → node 1
        face[1] = 2;  // face 2 → node 2
    }
    else if (dim == 2) {
        if (elemtype == 0) {  // triangle
            //face = (int*) malloc(sizeof(int) * (*nfe) * (nvf));
            // Each column corresponds to a face; column-major layout
            face[0 + 0*(nvf)] = 2;
            face[1 + 0*(nvf)] = 3;

            face[0 + 1*(nvf)] = 3;
            face[1 + 1*(nvf)] = 1;

            face[0 + 2*(nvf)] = 1;
            face[1 + 2*(nvf)] = 2;
        }
        else if (elemtype == 1) {  // quadrilateral
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(nvf)] = 1;
            face[1 + 0*(nvf)] = 2;

            face[0 + 1*(nvf)] = 2;
            face[1 + 1*(nvf)] = 3;

            face[0 + 2*(nvf)] = 3;
            face[1 + 2*(nvf)] = 4;

            face[0 + 3*(nvf)] = 4;
            face[1 + 3*(nvf)] = 1;
        }
    }
    else if (dim == 3) {
        if (elemtype == 0) {  // tetrahedron
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(nvf)] = 2;
            face[1 + 0*(nvf)] = 3;
            face[2 + 0*(nvf)] = 4;

            face[0 + 1*(nvf)] = 1;
            face[1 + 1*(nvf)] = 4;
            face[2 + 1*(nvf)] = 3;

            face[0 + 2*(nvf)] = 1;
            face[1 + 2*(nvf)] = 2;
            face[2 + 2*(nvf)] = 4;

            face[0 + 3*(nvf)] = 1;
            face[1 + 3*(nvf)] = 3;
            face[2 + 3*(nvf)] = 2;
        }
        else if (elemtype == 1) {  // hexahedron
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(nvf)] = 1;
            face[1 + 0*(nvf)] = 4;
            face[2 + 0*(nvf)] = 3;
            face[3 + 0*(nvf)] = 2;

            face[0 + 1*(nvf)] = 5;
            face[1 + 1*(nvf)] = 6;
            face[2 + 1*(nvf)] = 7;
            face[3 + 1*(nvf)] = 8;

            face[0 + 2*(nvf)] = 1;
            face[1 + 2*(nvf)] = 2;
            face[2 + 2*(nvf)] = 6;
            face[3 + 2*(nvf)] = 5;

            face[0 + 3*(nvf)] = 3;
            face[1 + 3*(nvf)] = 4;
            face[2 + 3*(nvf)] = 8;
            face[3 + 3*(nvf)] = 7;

            face[0 + 4*(nvf)] = 2;
            face[1 + 4*(nvf)] = 3;
            face[2 + 4*(nvf)] = 7;
            face[3 + 4*(nvf)] = 6;

            face[0 + 5*(nvf)] = 4;
            face[1 + 5*(nvf)] = 1;
            face[2 + 5*(nvf)] = 5;
            face[3 + 5*(nvf)] = 8;
        }
    }
    else {
        fprintf(stderr, "Error: Only dim = 1, 2, or 3 supported.\n");
        exit(1);
    }
            
    for (int i=0; i<nvf*nfe; i++) face[i] -= 1;
}

template<typename T>
void xiny2(int* out, const T* A, const T* B, int m, int n, int dim, double tol = 1e-12) {
    for (int i = 0; i < m; ++i) {
        out[i] = -1;
        for (int j = 0; j < n; ++j) {
            bool match = true;
            for (int d = 0; d < dim; ++d) {
                if (std::abs(A[i + n*d] - B[j + m*d]) > tol) {
                    match = false;
                    break;
                }
            }
            if (match) {
                out[i] = j;
                break;
            }
        }
    }
}

int permindex(int*& permind, const double* plocfc, int npf, int dim, int elemtype) 
{
    int ncols_out = 1;
    if (dim == 1) {         
        TemplateMalloc(&permind, 1, 0);
        permind[0] = 0;
    } 
    else if (dim == 2) {         
        TemplateMalloc(&permind, npf, 0);
        for (int i = 0; i < npf; ++i)
            permind[i] = npf - i - 1;
    } 
    else if (dim == 3 && elemtype == 0) {
        ncols_out = 3;        
        TemplateMalloc(&permind, npf*3, 0);
        double* plocfc2 = (double*) malloc(sizeof(double) * npf * 2);

        // [1 3 2]: swap columns
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 0*npf];
        }
        xiny2<double>(&permind[0], plocfc, plocfc2, npf, npf, 2, 1e-8);

        // [2 1 3]: 1 - xi - eta in col 1
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 0*npf] - plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 1*npf];
        }
        xiny2<double>(&permind[npf], plocfc, plocfc2, npf, npf, 2, 1e-8);

        // [3 2 1]: 1 - xi - eta in col 2
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 0*npf] - plocfc[i + 1*npf];
        }
        xiny2<double>(&permind[2*npf], plocfc, plocfc2, npf, npf, 2, 1e-8);

        CPUFREE(plocfc2);
    } 
    else if (dim == 3 && elemtype == 1) {
        ncols_out = 4;        
        TemplateMalloc(&permind, npf*4, 0);
        double* plocfc2 = (double*) malloc(sizeof(double) * npf * 2);

        // [1 4 3 2]: swap columns
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 0*npf];
        }
        xiny2<double>(&permind[0], plocfc, plocfc2, npf, npf, 2, 1e-8);
        
        // [2 1 4 3]: eta = 1 - eta
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 1*npf];
        }
        xiny2<double>(&permind[npf], plocfc, plocfc2, npf, npf, 2, 1e-8);

        // [3 2 1 4]: xi = 1 - eta, eta = 1 - xi
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 0*npf];
        }
        xiny2<double>(&permind[2*npf], plocfc, plocfc2, npf, npf, 2, 1e-8);

        // [4 3 2 1]: xi = 1 - xi
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = plocfc[i + 1*npf];
        }
        xiny2<double>(&permind[3*npf], plocfc, plocfc2, npf, npf, 2, 1e-8);

        CPUFREE(plocfc2);
    }    
        
    return ncols_out;
}

void compare(int a, int b, string s)
{
    cout<<s<<": ("<<a<<", " <<b<<")\n";
}

void buildConn(meshstruct& mesh, solstruct& sol, const appstruct& app, const masterstruct& master, 
               int* ti, int* boundaryConditions, int* intepartpts, int nbcm)
{                  
    int nd = master.ndims[0];     // spatial dimension    
    int elemtype = master.ndims[1]; 
    int npe = master.ndims[5]; // number of nodes on master element
    int npf = master.ndims[6]; // number of nodes on master face           
    int dim = nd;

    int nve = (elemtype==0) ? (dim + 1) : std::pow(2, dim);    
    int nvf = (dim == 3) ? (dim + elemtype) : dim;     
    int nfe = dim + (dim-1)*elemtype + 1;

    int ne = mesh.ndims[1];
    int nebmax = 4096; //mesh.ndims[6]; // maximum number of elements per block    
    int nfbmax = 8192; //mesh.ndims[8]; // maximum number of faces per block     

    int hybrid = app.problem[0];
    int mpiprocs = app.comm[1];
    int coupledinterface = app.problem[28];

    int ne1 = ne;    
    if (mpiprocs>1) {
      ne1 = mesh.elempartpts[0] + mesh.elempartpts[1];
      if (hybrid==0) ne1 += mesh.elempartpts[2];        
    }
        
    int* localfaces; 
    TemplateMalloc(&localfaces, nvf * nfe, 0);   
    getelemface(localfaces, dim, elemtype);

    int* permind;
    int ncols_out = permindex(permind, master.xpf, npf, dim, elemtype);

    int* facepartpts=nullptr;
    int* facepartbnd=nullptr;    
    int  sizes[10];
    int nf = connectivity(mesh.elemcon, mesh.facecon, mesh.f2e, facepartpts, 
              facepartbnd, sizes, mesh.bf, ti, mesh.elempart, localfaces, mesh.perm, 
              permind, boundaryConditions, dim, elemtype, nve, nvf, nfe, 
              npf, npe, ne, ne1, nbcm);       

    int nbe = 0;
    if (mpiprocs==1) {
      if (coupledinterface>0) {        
        vector<int> me(3), tm(2); 
        me[0] = 0; me[1] = intepartpts[0]; me[2] = me[1] + intepartpts[1];
        tm[0] = 0; tm[1] = -1;
        nbe = mkfaceblocks(mesh.eblks, me.data(), tm.data(), 3, nebmax); 
      } else {
        nbe = divide_interval(mesh.eblks, ne, nebmax);
      }
    } else {
      if (coupledinterface>0) {        
        vector<int> me(5), tm(4); 
        me[0] = 0; me[1] = intepartpts[0]; me[2] = me[1] + intepartpts[1]; me[3] = me[2] + intepartpts[2]; me[4] = me[3] + intepartpts[3];
        tm[0] = 0; tm[1] = -1; tm[2] = 1; tm[3] = 2;
        nbe = mkfaceblocks(mesh.eblks, me.data(), tm.data(), 5, nebmax); 
      } else {
        vector<int> me(4), tm(3); 
        me[0] = 0; me[1] = mesh.elempartpts[0]; me[2] = me[1] + mesh.elempartpts[1]; me[3] = me[2] + mesh.elempartpts[2];
        tm[0] = 0; tm[1] = 1; tm[2] = 2;
        nbe = mkfaceblocks(mesh.eblks, me.data(), tm.data(), 4, nebmax); 
      }      
    }

    int nbc = sizes[3];
    int n = 1 + nbc + 1;
    vector<int> mf(n); mf[0] = 0; 
    for(int i=1; i<n; i++) mf[i] = mf[i-1] + facepartpts[i-1];
    int nbf = mkfaceblocks(mesh.fblks, mf.data(), facepartbnd, n, nfbmax); 

    int neb = 0, nfb = 0;
    for (int i=0; i<nbe; i++) {
      int k = mesh.eblks[3*i + 1] - mesh.eblks[3*i + 0] + 1;
      neb = (k > neb) ? k : neb;
    }
    for (int i=0; i<nbf; i++) {
      int k = mesh.fblks[3*i + 1] - mesh.fblks[3*i + 0] + 1;
      nfb = (k > nfb) ? k : nfb;
    }

    int nfc = npf * nf;
    vector<int> facecon1(nfc);
    vector<int> facecon2(nfc);
    for (int i=0; i<nfc; i++) {
      facecon1[i] = mesh.facecon[0 + 2*i];
      facecon2[i] = mesh.facecon[1 + 2*i];
    }      

    int new_nf = removeBoundaryFaces(facecon2.data(), mesh.fblks, nbf, npf, nf);     
    int nent1 = mkdge2dgf(mesh.rowe2f1, mesh.cole2f1, mesh.ent2ind1, facecon1.data(), npf*nf, npe*ne); 
    int nent2 = mkdge2dgf(mesh.rowe2f2, mesh.cole2f2, mesh.ent2ind2, facecon2.data(), npf*new_nf, npe*ne); 

    dstype *xcg=nullptr;
    TemplateMalloc(&xcg, npe*dim*ne, 0);
    TemplateMalloc(&mesh.cgelcon, npe*ne, 0);
    int ncgnodes = mkelconcg_hashgrid(xcg, mesh.cgelcon, sol.xdg, npe, dim, ne);
    int ncgdof = mkent2elem(mesh.rowent2elem, mesh.colent2elem, mesh.cgelcon, npe, ne);
    map_cgent2dgent(mesh.cgent2dgent, mesh.rowent2elem, mesh.colent2elem, xcg, sol.xdg,  npe, dim, ncgdof);                                          
    TemplateMalloc(&sol.xcg, dim*ncgnodes, 0);
    for (int i=0; i<dim*ncgnodes; i++) sol.xcg[i] = xcg[i];    
    sol.szxcg = dim*ncgnodes;

    // compare(mesh.ndims[0], dim, "dim");
    // compare(mesh.ndims[1], ne, "ne");
    // compare(mesh.ndims[2], nf, "nf");             // number of faces in this subdomain
    // //compare(mesh.ndims[3], ncgnodes, "ncgnodes"); // number of vertices in this subdomain
    // compare(mesh.ndims[4], nfe, "nfe");           // number of faces per element
    // compare(mesh.ndims[5], nbe, "nbe");           // number of blocks for elements
    // compare(mesh.ndims[6], neb, "neb");           // max number of elements per block
    // compare(mesh.ndims[7], nbf, "nbf");           // number of blocks for faces
    // compare(mesh.ndims[8], nfb, "nfb");           // max number of faces per block
    // 
    // compare(mesh.nsize[1], 2 * npf * nf, "mesh.nsize[1]");
    // compare(mesh.nsize[2], 3 * nbe, "mesh.nsize[2]");
    // compare(mesh.nsize[3], 3 * nbf, "mesh.nsize[3]");
    // compare(mesh.nsize[11], npe * ne, "mesh.nsize[11]");
    // compare(mesh.nsize[12], ncgdof + 1, "mesh.nsize[12]");
    // compare(mesh.nsize[13], mesh.rowent2elem[ncgdof], "mesh.nsize[13]");
    // compare(mesh.nsize[14], mesh.rowent2elem[ncgdof], "mesh.nsize[14]");
    // compare(mesh.nsize[15], nent1 + 1, "mesh.nsize[15]");
    // compare(mesh.nsize[16], mesh.rowe2f1[nent1], "mesh.nsize[16]");
    // compare(mesh.nsize[17], npe * ne, "mesh.nsize[17]");
    // compare(mesh.nsize[18], nent2 + 1, "mesh.nsize[18]");
    // compare(mesh.nsize[19], mesh.rowe2f2[nent2], "mesh.nsize[19]");
    // compare(mesh.nsize[20], npe * ne, "mesh.nsize[20]");
    // compare(mesh.nsize[21], 4 * nf, "mesh.nsize[21]");
    
    mesh.ndims[0] = dim;
    mesh.ndims[1] = ne;
    mesh.ndims[2] = nf; // number of faces in this subdomain 
    //mesh.ndims[3] = ncgnodes; // number of vertices in this subdomain       
    mesh.ndims[4] = nfe; // number of faces per element            
    mesh.ndims[5] = nbe; // number of blocks for elements 
    mesh.ndims[6] = neb; // maximum number of elements per block
    mesh.ndims[7] = nbf; // number of blocks for faces   
    mesh.ndims[8] = nfb; // maximum number of faces per block     
    
    mesh.szfacecon = mesh.nsize[1] = 2 * npf * nf;
    mesh.szeblks = mesh.nsize[2] = 3 * nbe;
    mesh.szfblks = mesh.nsize[3] = 3 * nbf;
    mesh.szcgelcon = mesh.nsize[11] = npe * ne;
    mesh.szrowent2elem = mesh.nsize[12] = ncgdof + 1;
    mesh.szcgent2dgent = mesh.nsize[13] = mesh.rowent2elem[ncgdof];
    mesh.szcolent2elem = mesh.nsize[14] = mesh.rowent2elem[ncgdof];
    mesh.szrowe2f1 = mesh.nsize[15] = nent1 + 1;
    mesh.szcole2f1 = mesh.nsize[16] = mesh.rowe2f1[nent1];
    mesh.szent2ind1 = mesh.nsize[17] = npe*ne;
    mesh.szrowe2f2 = mesh.nsize[18] = nent2 + 1;
    mesh.szcole2f2 = mesh.nsize[19] = mesh.rowe2f2[nent2];
    mesh.szent2ind2 = mesh.nsize[20] = npe*ne;
    mesh.szf2e = mesh.nsize[21] = 4*nf;    
    mesh.szelemcon = mesh.nsize[22] = npf * nfe * ne;
    
    CPUFREE(xcg);
    CPUFREE(localfaces);
    CPUFREE(permind);
    CPUFREE(facepartpts);
    CPUFREE(facepartbnd);
}

struct Conn 
{    
    dstype *xcg;
    int *elemcon, *facecon, *f2e, *eblks, *fblks;
    int *cgelcon, *rowent2elem, *colent2elem, *cgent2dgent;
    int *rowe2f1, *cole2f1, *ent2ind1, *rowe2f2, *cole2f2, *ent2ind2;
    //int ncgnodes, ncgdof, nbe, nbf, neb, nfb, nf;
};

void buildConn(Conn& conn, meshstruct& mesh, solstruct& sol, const appstruct& app, const masterstruct& master, 
               int* ti, int* boundaryConditions, int* intepartpts, int nbcm)
{                  
    int nd = master.ndims[0];     // spatial dimension    
    int elemtype = master.ndims[1]; 
    int npe = master.ndims[5]; // number of nodes on master element
    int npf = master.ndims[6]; // number of nodes on master face           
    int dim = nd;

    int nve = (elemtype==0) ? (dim + 1) : std::pow(2, dim);    
    int nvf = (dim == 3) ? (dim + elemtype) : dim;     
    int nfe = dim + (dim-1)*elemtype + 1;

    int ne = mesh.ndims[1];
    int nebmax = mesh.ndims[6]; // maximum number of elements per block    
    int nfbmax = mesh.ndims[8]; // maximum number of faces per block     

    int hybrid = app.problem[0];
    int mpiprocs = app.comm[1];
    int coupledinterface = app.problem[28];

    int ne1 = ne;    
    if (mpiprocs>1) {
      ne1 = mesh.elempartpts[0] + mesh.elempartpts[1];
      if (hybrid==0) ne1 += mesh.elempartpts[2];        
    }
        
    int* localfaces; 
    TemplateMalloc(&localfaces, nvf * nfe, 0);   
    getelemface(localfaces, dim, elemtype);

    int* permind;
    int ncols_out = permindex(permind, master.xpf, npf, dim, elemtype);

    int* facepartpts=nullptr;
    int* facepartbnd=nullptr;    
    int  sizes[10];
    int nf = connectivity(conn.elemcon, conn.facecon, conn.f2e, facepartpts, 
              facepartbnd, sizes, mesh.bf, ti, mesh.elempart, localfaces, mesh.perm, 
              permind, boundaryConditions, dim, elemtype, nve, nvf, nfe, 
              npf, npe, ne, ne1, nbcm);       

    //print2iarray(conn.f2e, nfe, nf);

    int nbe = 0;
    if (mpiprocs==1) {
      if (coupledinterface>0) {        
        vector<int> me(3), tm(2); 
        me[0] = 0; me[1] = intepartpts[0]; me[2] = me[1] + intepartpts[1];
        tm[0] = 0; tm[1] = -1;
        nbe = mkfaceblocks(conn.eblks, me.data(), tm.data(), 3, nebmax); 
      } else {
        nbe = divide_interval(conn.eblks, ne, nebmax);
      }
    } else {
      if (coupledinterface>0) {        
        vector<int> me(5), tm(4); 
        me[0] = 0; me[1] = intepartpts[0]; me[2] = me[1] + intepartpts[1]; me[3] = me[2] + intepartpts[2]; me[4] = me[3] + intepartpts[3];
        tm[0] = 0; tm[1] = -1; tm[2] = 1; tm[3] = 2;
        nbe = mkfaceblocks(conn.eblks, me.data(), tm.data(), 5, nebmax); 
      } else {
        vector<int> me(4), tm(3); 
        me[0] = 0; me[1] = mesh.elempartpts[0]; me[2] = me[1] + mesh.elempartpts[1]; me[3] = me[2] + mesh.elempartpts[2];
        tm[0] = 0; tm[1] = 1; tm[2] = 2;
        nbe = mkfaceblocks(conn.eblks, me.data(), tm.data(), 4, nebmax); 
      }      
    }

    int nbc = sizes[3];
    int n = 1 + nbc + 1;
    vector<int> mf(n); mf[0] = 0; 
    for(int i=1; i<n; i++) mf[i] = mf[i-1] + facepartpts[i-1];
    int nbf = mkfaceblocks(conn.fblks, mf.data(), facepartbnd, n, nfbmax); 

    int neb = 0, nfb = 0;
    for (int i=0; i<nbe; i++) {
      int k = conn.eblks[3*i + 1] - conn.eblks[3*i + 0] + 1;
      neb = (k > neb) ? k : neb;
    }
    for (int i=0; i<nbf; i++) {
      int k = conn.fblks[3*i + 1] - conn.fblks[3*i + 0] + 1;
      nfb = (k > nfb) ? k : nfb;
    }

    int nfc = npf * nf;
    vector<int> facecon1(nfc);
    vector<int> facecon2(nfc);
    for (int i=0; i<nfc; i++) {
      facecon1[i] = mesh.facecon[0 + 2*i];
      facecon2[i] = mesh.facecon[1 + 2*i];
    }      

    int new_nf = removeBoundaryFaces(facecon2.data(), conn.fblks, nbf, npf, nf);     
    int nent1 = mkdge2dgf(conn.rowe2f1, conn.cole2f1, conn.ent2ind1, facecon1.data(), npf*nf, npe*ne); 
    int nent2 = mkdge2dgf(conn.rowe2f2, conn.cole2f2, conn.ent2ind2, facecon2.data(), npf*new_nf, npe*ne); 

    dstype *xcg=nullptr;
    TemplateMalloc(&xcg, npe*dim*ne, 0);
    TemplateMalloc(&conn.cgelcon, npe*ne, 0);
    int ncgnodes = mkelconcg_hashgrid(xcg, conn.cgelcon, sol.xdg, npe, dim, ne);
    int ncgdof = mkent2elem(conn.rowent2elem, conn.colent2elem, conn.cgelcon, npe, ne);
    map_cgent2dgent(conn.cgent2dgent, conn.rowent2elem, conn.colent2elem, xcg, sol.xdg,  npe, dim, ncgdof);                                          
    TemplateMalloc(&conn.xcg, dim*ncgnodes, 0);
    for (int i=0; i<dim*ncgnodes; i++) conn.xcg[i] = xcg[i];
    
    compare(mesh.ndims[0], dim, "dim");
    compare(mesh.ndims[1], ne, "ne");
    compare(mesh.ndims[2], nf, "nf");             // number of faces in this subdomain
    //compare(mesh.ndims[3], ncgnodes, "ncgnodes"); // number of vertices in this subdomain
    compare(mesh.ndims[4], nfe, "nfe");           // number of faces per element
    compare(mesh.ndims[5], nbe, "nbe");           // number of blocks for elements
    compare(mesh.ndims[6], neb, "neb");           // max number of elements per block
    compare(mesh.ndims[7], nbf, "nbf");           // number of blocks for faces
    compare(mesh.ndims[8], nfb, "nfb");           // max number of faces per block
       
    compare(mesh.nsize[1], 2 * npf * nf, "mesh.nsize[1]");
    compare(mesh.nsize[2], 3 * nbe, "mesh.nsize[2]");
    compare(mesh.nsize[3], 3 * nbf, "mesh.nsize[3]");
    compare(mesh.nsize[11], npe * ne, "mesh.nsize[11]");
    compare(mesh.nsize[12], ncgdof + 1, "mesh.nsize[12]");
    compare(mesh.nsize[13], mesh.rowent2elem[ncgdof], "mesh.nsize[13]");
    compare(mesh.nsize[14], mesh.rowent2elem[ncgdof], "mesh.nsize[14]");
    compare(mesh.nsize[15], nent1 + 1, "mesh.nsize[15]");
    compare(mesh.nsize[16], mesh.rowe2f1[nent1], "mesh.nsize[16]");
    compare(mesh.nsize[17], npe * ne, "mesh.nsize[17]");
    compare(mesh.nsize[18], nent2 + 1, "mesh.nsize[18]");
    compare(mesh.nsize[19], mesh.rowe2f2[nent2], "mesh.nsize[19]");
    compare(mesh.nsize[20], npe * ne, "mesh.nsize[20]");
    compare(mesh.nsize[21], 4 * nf, "mesh.nsize[21]");
    compare(mesh.nsize[22], npf * nfe * ne, "mesh.nsize[21]");
    
    CPUFREE(xcg);
    CPUFREE(localfaces);
    CPUFREE(permind);
    CPUFREE(facepartpts);
    CPUFREE(facepartbnd);
}

void maxdiff(const dstype *a, const dstype *b, int n, string s)
{
    dstype e = 0.0;
    for (int i = 0; i < n; i++) {
        dstype d = fabs(a[i] - b[i]);
        e = (e > d) ? e : d;
    }
    cout<<"Maximum error in "<<s<<" = "<<e<<", length = "<<n<<endl;
}

void maxdiff(const int *a, const int *b, int n, string s)
{
    int e = 0;
    for (int i = 0; i < n; i++) {
        int d = abs(a[i] - b[i]);
        e = (e > d) ? e : d;
    }
    cout<<"Maximum error in "<<s<<" = "<<e<<", length = "<<n<<endl;
}

void checkConn(meshstruct& mesh, solstruct& sol, const appstruct& app, const masterstruct& master, 
               int* ti, int* boundaryConditions, int* intepartpts, int nbcm)
{
    Conn conn;
    buildConn(conn, mesh, sol, app, master, ti, boundaryConditions, intepartpts, nbcm);

    maxdiff(conn.xcg, sol.xcg, sol.szxcg, "xcg");

    maxdiff(conn.elemcon, mesh.elemcon, mesh.szelemcon, "elemcon");
    maxdiff(conn.facecon, mesh.facecon, mesh.szfacecon, "facecon");    
    maxdiff(conn.f2e,     mesh.f2e,     mesh.szf2e,     "f2e");

    maxdiff(conn.eblks,   mesh.eblks,   mesh.szeblks,   "eblks");
    maxdiff(conn.fblks,   mesh.fblks,   mesh.szfblks,   "fblks");
    
    maxdiff(conn.cgelcon,    mesh.cgelcon,    mesh.szcgelcon,    "cgelcon");
    maxdiff(conn.rowent2elem, mesh.rowent2elem, mesh.szrowent2elem, "rowent2elem");
    maxdiff(conn.colent2elem, mesh.colent2elem, mesh.szcolent2elem, "colent2elem");
    maxdiff(conn.cgent2dgent, mesh.cgent2dgent, mesh.szcgent2dgent, "cgent2dgent");
    
    maxdiff(conn.rowe2f1,  mesh.rowe2f1,  mesh.szrowe2f1,  "rowe2f1");
    maxdiff(conn.cole2f1,  mesh.cole2f1,  mesh.szcole2f1,  "cole2f1");
    maxdiff(conn.ent2ind1, mesh.ent2ind1, mesh.szent2ind1, "ent2ind1");
    
    maxdiff(conn.rowe2f2,  mesh.rowe2f2,  mesh.szrowe2f2,  "rowe2f2");
    maxdiff(conn.cole2f2,  mesh.cole2f2,  mesh.szcole2f2,  "cole2f2");
    maxdiff(conn.ent2ind2, mesh.ent2ind2, mesh.szent2ind2, "ent2ind2");    
}

#endif

