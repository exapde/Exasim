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

struct Conn 
{    
    vector<double> cgnodes;
    vector<int> elemcon, facecon, bf, f2t, facepartpts, facepartbnd, eblks, fblks;
    vector<int> cgelcon, rowent2elem, colent2elem, cgent2dgent;
    vector<int> rowe2f1, cole2f1, ent2ind1, rowe2f2, cole2f2, ent2ind2;
    int ncgnodes=0, ncgdof=0, nbe=0, nbf=0, neb=0, nfb=0, nf=0;
};
    
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
        //print2iarray(perm, npf, nfe);
        //print2iarray(permind, npf, nvf);
        
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
                    print2iarray(f1, 1, nvf);
                    print2iarray(f2, 1, nvf);
                    print2iarray(face, nvf, nfe);
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

int connectivity(vector<int>& elemcon, vector<int>& facecon, vector<int>& bf, vector<int>& f2t, 
                 vector<int>& facepartpts, vector<int>& facepartbnd, 
                 const int* t, const int* f, const int* elempart, const int* localfaces,  
                 const int* perm, const int* permind,  const int* bcm, int dim, int elemtype, 
                 int nve, int nvf, int nfe, int npf, int npe, int ne, int ne1, int nbcm) 
{      
      bf.resize(nfe * ne);
  
      int* fi = (int*)malloc(nfe * ne * sizeof(int));      
      select_columns(fi, f, elempart, nfe, ne);       
      apply_bcm(bf.data(), fi, bcm, nfe*ne, nbcm);                   
      
//       print2iarray(fi, nfe, ne);
//       print2iarray(bf.data(), nfe, ne);
//       print2iarray(bcm, 1, nbcm);
//       error("here");
      
      CPUFREE(fi);
                                                      
      int* ti = (int*)malloc(nve * ne * sizeof(int));                  
      select_columns(ti, t, elempart, nve, ne); 
      
      int* f2t_tmp = (int*)malloc(4 * nfe * ne * sizeof(int));
      int nf1 = mkf2e_hash(f2t_tmp, ti, localfaces, ne, nve, nvf, nfe);       
      
//       print2iarray(elempart, 1, ne);
//       print2iarray(ti, nve, ne);
//       print2iarray(f2t_tmp, 4, nf1);
      
      int* ind = (int*)malloc(nf1 * sizeof(int));
      int nf = 0;
      for (int i=0; i<nf1; i++)
        if (f2t_tmp[0 + 4*i] < ne1)
          ind[nf++] = i;                
      
      f2t.resize(4*nf);
      select_columns(f2t.data(), f2t_tmp, ind, 4, nf);        
      CPUFREE(f2t_tmp);
      
//       print2iarray(f2t.data(), 4, nf);
//       error("here");
      
      int nifaces = find(f2t.data(), 0, 4, nf, 2, 2);
      int nbfaces = find(f2t.data(), -1, 4, nf, 2, 0);      
      int* ifaces = (int*)malloc(nifaces * sizeof(int));
      int* bfaces = (int*)malloc(nbfaces * sizeof(int));
      find(ifaces, f2t.data(), 0, 4, nf, 2, 2);  // interior faces
      find(bfaces, f2t.data(), -1, 4, nf, 2, 0); // boundary faces
            
      //print2iarray(ifaces, 1, nifaces);
      //print2iarray(bfaces, 1, nbfaces);
      //cout<<nifaces<<", "<<nbfaces<<endl;
      //int find(int* indices, const int* a, int b, int m, int n, int k, int opts) 
      
      int* bfa = (int*)malloc(nbfaces * sizeof(int));
      for (int i=0; i<nbfaces; i++) {
        int e = f2t[0 + 4*bfaces[i]];
        int l = f2t[1 + 4*bfaces[i]];
        bfa[i] = bf[l + nfe*e];
      }                    
                  
      int* sortedbfaces = (int*)malloc(nbfaces * sizeof(int));
      simple_bubble_sort(sortedbfaces, ind, bfa, nbfaces);                  
      
//       print2iarray(bfa, 1, nbfaces);
//       print2iarray(ind, 1, nbfaces);
      
      facepartbnd.resize(nbcm+1);
      facepartpts.resize(nbcm+1);
      facepartbnd[0] = 0;
      facepartpts[0] = nifaces;
      int nbc = unique_count(&facepartbnd[1], &facepartpts[1], sortedbfaces, nbfaces);      
      facepartpts.erase(facepartpts.begin() + 1 + nbc, facepartpts.end());
      facepartbnd.erase(facepartbnd.begin() + 1 + nbc, facepartbnd.end());
      
//       print2iarray(facepartpts.data(), 1, facepartpts.size());
//       print2iarray(facepartbnd.data(), 1, facepartbnd.size());
      
      int* bind = (int*)malloc(nf * sizeof(int));
      for (int i=0; i<nifaces; i++) bind[i] = ifaces[i];
      for (int i=0; i<nbfaces; i++) bind[nifaces+i] = bfaces[ind[i]];
      permute_columns(f2t.data(), bind, 4, nf); 
      
      CPUFREE(ifaces); CPUFREE(bfaces); CPUFREE(bfa); CPUFREE(ind); CPUFREE(bind);  CPUFREE(sortedbfaces);                
      
      //print2iarray(bfaces, 1, nbfaces);
      //print2iarray(ind, 1, nbfaces);
      //print2iarray(bind, 1, nf);
      //print2iarray(f2t.data(), 4, nf);
      
      fix_f2t_consistency(f2t.data(), elempart, nf); 
      
      //print2iarray(f2t.data(), 4, nf);
      
      elemcon.resize(npf * nfe * ne);      
      std::fill(elemcon.begin(), elemcon.end(), -1);
      facecon.resize(npf * 2 * nf);      
      std::fill(facecon.begin(), facecon.end(), -1);
      build_connectivity(elemcon.data(), facecon.data(), f2t.data(), localfaces, ti, perm, 
              permind, dim, elemtype, npf, nfe, npe, ne, nf); 
      
      //print2iarray(elemcon.data(), npf, nfe);            
      
//       print2iarray(perm, 1, npf);
//       print2iarray(permind, 1, npf);
      
      CPUFREE(ti);
      
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

// Converts element-to-node connectivity (cgelcon) to node-to-element connectivity (CRS format)
// cgelcon: [nrow x ne] column-major, zero-based node indices
int mkent2elem(
    vector<int>& rowent2elem, // output [ndof+1], allocated inside
    vector<int>& colent2elem,  // output [nnz], allocated inside        
    vector<int>& cgelcon, // [nrow * ne], column-major
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
    rowent2elem.resize(ndof+1);
    std::fill(rowent2elem.begin(), rowent2elem.end(), 0);

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

    colent2elem.resize(rowent2elem[ndof]);
    std::fill(colent2elem.begin(), colent2elem.end(), 0);
    
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
    vector<int>& cgent2dgent,        // [nnz], output mapping (same shape as colent2elem)
    const vector<int>& rowent2elem, // [nent+1]
    const vector<int>& colent2elem, // [nnz]    
    const vector<double>& cgnodes,  // [nent * dim], row-major
    const vector<double>& dgnodes,  // [npe * dim * ne], column-major
    int npe, int dim, int nent) 
{    
    cgent2dgent.resize(rowent2elem[nent]);
    for (int i = 0; i < nent; ++i) {
        int start = rowent2elem[i];
        int end = rowent2elem[i + 1];
        const double* xcg = &cgnodes[i * dim];
        for (int j = start; j < end; ++j) {
            int elem = colent2elem[j];
            const double* xdg = &dgnodes[npe * dim * elem];
            int in = xiny(xcg, xdg, npe, dim);
            if (in==-1) {
              print2darray(xcg, 1, dim);
              print2darray(xdg, npe, dim);
              for (int k = 0; k < npe; ++k) 
                for (int d = 0; d < dim; ++d) 
                  printf("%d %d: %g\n", k, d, fabs(xcg[d] - xdg[k + npe * d]));
              error("xcg and xdg do not match.");            
            }        
            cgent2dgent[j] = elem * npe + in; // global DG index
        }
    }
}

// Construct CRS-style mapping from facecon to per-entity degrees of freedom
int mkdge2dgf(vector<int>& rowdge2dgf, vector<int>& coldge2dgf, vector<int>& ent2ind,
              const vector<int>& facecon, int ndgf, int entmax) 
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
    ent2ind.resize(entmax);
    std::fill(ent2ind.begin(), ent2ind.end(), -1);
    for (int i = 0; i < nent; ++i)
        ent2ind[tmp[i]] = i;

    // Step 4: Count occurrences per entity
    //int* rowdge2dgf = (int*)calloc(nent + 1, sizeof(int));
    rowdge2dgf.resize(nent+1);
    std::fill(rowdge2dgf.begin(), rowdge2dgf.end(), 0);
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
    coldge2dgf.resize(total);
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

int removeBoundaryFaces(std::vector<int>& facecon, const std::vector<int>& fblks, int m, int npf, int nf) 
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
    facecon = std::move(new_facecon);
    
    return new_nf;
}

int divide_interval(vector<int>& intervals, int n, int m) 
{
    if (n <= 0 || m <= 0) return 0;

    int num_intervals = (n + m - 1) / m; // ceil(n/m)
    intervals.resize(3*num_intervals);
            
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

int mkfaceblocks(vector<int>& nm, const vector<int>& mf, const vector<int>& bcm, int nmf_len, int ns) 
{
    if (ns <= 0) ns = 2048;  // default value

    int max_blocks = 0;
    for (int i = 0; i < nmf_len - 1; ++i)
        max_blocks += (mf[i + 1] - mf[i] + ns - 1) / ns;

    nm.resize(3 * max_blocks);
    int count = 0;

    for (int i = 0; i < nmf_len - 1; ++i) {
        int nf = mf[i + 1] - mf[i];
        vector<int> intervals; 
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

void buildConn(Conn& conn, const PDE& pde, const Mesh& mesh, const Master& master, const DMD& dmd)
{                  
    //int ne = mesh.ne;
    //int ne1 = mesh.ne;   
    int ne = dmd.elempart.size(); 
    int ne1 = ne;
    if (pde.mpiprocs>1) {      
      ne1 = dmd.elempartpts[0] + dmd.elempartpts[1];
      if (pde.hybrid==0) ne1 += dmd.elempartpts[2];        
    }
    
    conn.nf = connectivity(conn.elemcon, conn.facecon, conn.bf, conn.f2t, conn.facepartpts, conn.facepartbnd, 
              mesh.t.data(), mesh.f.data(), dmd.elempart.data(), mesh.localfaces.data(), master.perm.data(), 
              master.permind.data(), mesh.boundaryConditions.data(), mesh.dim, mesh.elemtype, mesh.nve, 
              mesh.nvf, mesh.nfe, master.npf, master.npe, ne, ne1, mesh.nbcm);       
    
    cout << "Finished connectivity" << endl;

    //writesol(pde, mesh, master, dmd.elempart, conn.facepartpts, 0, pde.datapath + "/sol.bin");        
    //print2iarray(conn.elemcon.data(), master.npf*mesh.nfe, ne);

    if (pde.mpiprocs==1) {
      if (pde.coupledinterface>0) {        
        vector<int> me(3), tm(2); 
        me[0] = 0; me[1] = dmd.intepartpts[0]; me[2] = me[1] + dmd.intepartpts[1];
        tm[0] = 0; tm[1] = -1;
        conn.nbe = mkfaceblocks(conn.eblks, me, tm, 3, pde.neb); 
      } else {
        conn.nbe = divide_interval(conn.eblks, ne, pde.neb);
      }
    } else {
      if (pde.coupledinterface>0) {        
        vector<int> me(5), tm(4); 
        me[0] = 0; me[1] = dmd.intepartpts[0]; me[2] = me[1] + dmd.intepartpts[1]; me[3] = me[2] + dmd.intepartpts[2]; me[4] = me[3] + dmd.intepartpts[3];
        tm[0] = 0; tm[1] = -1; tm[2] = 1; tm[3] = 2;
        conn.nbe = mkfaceblocks(conn.eblks, me, tm, 5, pde.neb); 
      } else {
        vector<int> me(4), tm(3); 
        me[0] = 0; me[1] = dmd.elempartpts[0]; me[2] = me[1] + dmd.elempartpts[1]; me[3] = me[2] + dmd.elempartpts[2];
        tm[0] = 0; tm[1] = 1; tm[2] = 2;
        conn.nbe = mkfaceblocks(conn.eblks, me, tm, 4, pde.neb); 
      }      
    }

    int n = 1 + conn.facepartpts.size();
    vector<int> mf(n); mf[0] = 0; 
    for(int i=1; i<n; i++) mf[i] = mf[i-1] + conn.facepartpts[i-1];
    conn.nbf = mkfaceblocks(conn.fblks, mf, conn.facepartbnd, n, pde.nfb); 
    
//     print2iarray(conn.facepartpts.data(), 1, conn.facepartpts.size());
//     print2iarray(conn.eblks.data(), 3, conn.nbe);
//     print2iarray(conn.fblks.data(), 3, conn.nbf);   
            
    conn.neb = 0, conn.nfb = 0;
    for (int i=0; i<conn.nbe; i++) {
      int k = conn.eblks[3*i + 1] - conn.eblks[3*i + 0] + 1;
      conn.neb = (k > conn.neb) ? k : conn.neb;
    }
    for (int i=0; i<conn.nbf; i++) {
      int k = conn.fblks[3*i + 1] - conn.fblks[3*i + 0] + 1;
      conn.nfb = (k > conn.nfb) ? k : conn.nfb;
    }

    //print2iarray(conn.eblks.data(), 3, nbe);
    //print2iarray(conn.fblks.data(), 3, nbf);
    
    int nfc = conn.facecon.size()/2;
    vector<int> facecon1(nfc);
    vector<int> facecon2(nfc);
    for (int i=0; i<nfc; i++) {
      facecon1[i] = conn.facecon[0 + 2*i];
      facecon2[i] = conn.facecon[1 + 2*i];
    }      
    
    int new_nf = removeBoundaryFaces(facecon2, conn.fblks, conn.nbf, master.npf, nfc/master.npf);     
    
    //print2iarray(facecon1.data(), master.npf, nfc/master.npf);
    //print2iarray(facecon2.data(), master.npf, new_nf);
    
    mkdge2dgf(conn.rowe2f1, conn.cole2f1, conn.ent2ind1, facecon1, facecon1.size(), master.npe*ne); 
    mkdge2dgf(conn.rowe2f2, conn.cole2f2, conn.ent2ind2, facecon2, facecon2.size(), master.npe*ne); 

    cout << "Finished mkdge2dgf" << endl;

    conn.cgnodes.resize(master.npe * mesh.dim * ne);
    conn.cgelcon.resize(master.npe * ne);
            
    if (pde.mpiprocs==1) {    
      conn.ncgnodes = mkelconcg_hashgrid(conn.cgnodes.data(), conn.cgelcon.data(), mesh.xdg.data(), master.npe, mesh.dim, ne);
      conn.ncgdof = mkent2elem(conn.rowent2elem, conn.colent2elem, conn.cgelcon, master.npe, ne);
      map_cgent2dgent(conn.cgent2dgent, conn.rowent2elem, conn.colent2elem, conn.cgnodes, mesh.xdg,  master.npe, mesh.dim, conn.ncgdof);                                      
    } else {
      vector<double> tm (master.npe*mesh.dim*dmd.elempart.size());
      select_columns(tm.data(), mesh.xdg.data(), dmd.elempart.data(), master.npe*mesh.dim, dmd.elempart.size());
      conn.ncgnodes = mkelconcg_hashgrid(conn.cgnodes.data(), conn.cgelcon.data(), tm.data(), master.npe, mesh.dim, ne);
      conn.ncgdof = mkent2elem(conn.rowent2elem, conn.colent2elem, conn.cgelcon, master.npe, ne);
      map_cgent2dgent(conn.cgent2dgent, conn.rowent2elem, conn.colent2elem, conn.cgnodes, tm,  master.npe, mesh.dim, conn.ncgdof);                                      
    }    
    conn.cgnodes.resize(mesh.dim * conn.ncgnodes);

    cout << "Finished map_cgent2dgent" << endl;

//     if (pde.debugmode==1 && pde.mpiprocs==1) {
//       //mesh.t.data(), mesh.f.data(), dmd.elempart.data(),
//       writearray2file(pde.datapath + "/t.bin", mesh.t.data(), mesh.t.size());
//       writearray2file(pde.datapath + "/f.bin", mesh.f.data(), mesh.f.size());
//       writearray2file(pde.datapath + "/elempart.bin", dmd.elempart.data(), dmd.elempart.size());
//       writearray2file(pde.datapath + "/localfaces.bin", mesh.localfaces.data(), mesh.localfaces.size());
//       writearray2file(pde.datapath + "/boundaryConditions.bin", mesh.boundaryConditions.data(), mesh.boundaryConditions.size());
//       writearray2file(pde.datapath + "/permind.bin", master.permind.data(), master.permind.size());              
//       writearray2file(pde.datapath + "/elemcon.bin", conn.elemcon.data(), conn.elemcon.size());
//       writearray2file(pde.datapath + "/facecon.bin", conn.facecon.data(), conn.facecon.size());
//       writearray2file(pde.datapath + "/bf.bin", conn.bf.data(), conn.bf.size());      
//       writearray2file(pde.datapath + "/f2t.bin", conn.f2t.data(), conn.f2t.size());      
//       writearray2file(pde.datapath + "/facepartpts.bin", conn.facepartpts.data(), conn.facepartpts.size());      
//       writearray2file(pde.datapath + "/facepartbnd.bin", conn.facepartbnd.data(), conn.facepartbnd.size());      
//       writearray2file(pde.datapath + "/eblks.bin", conn.eblks.data(), conn.eblks.size());      
//       writearray2file(pde.datapath + "/fblks.bin", conn.fblks.data(), conn.fblks.size());      
//       writearray2file(pde.datapath + "/facecon1.bin", facecon1.data(), facecon1.size());
//       writearray2file(pde.datapath + "/facecon2.bin", facecon2.data(), facecon2.size());
//       writearray2file(pde.datapath + "/rowe2f1.bin", conn.rowe2f1.data(), conn.rowe2f1.size());
//       writearray2file(pde.datapath + "/rowe2f2.bin", conn.rowe2f2.data(), conn.rowe2f2.size());
//       writearray2file(pde.datapath + "/cole2f1.bin", conn.cole2f1.data(), conn.cole2f1.size());
//       writearray2file(pde.datapath + "/cole2f2.bin", conn.cole2f2.data(), conn.cole2f2.size());
//       writearray2file(pde.datapath + "/ent2ind1.bin", conn.ent2ind1.data(), conn.ent2ind1.size());
//       writearray2file(pde.datapath + "/ent2ind2.bin", conn.ent2ind2.data(), conn.ent2ind2.size());
//       //writearray2file(pde.datapath + "/cgnodes.bin", conn.cgnodes.data(), conn.cgnodes.size());
//       writearray2file(pde.datapath + "/cgelcon.bin", conn.cgelcon.data(), conn.cgelcon.size());
//       writearray2file(pde.datapath + "/rowent2elem.bin", conn.rowent2elem.data(), conn.rowent2elem.size());
//       writearray2file(pde.datapath + "/colent2elem.bin", conn.colent2elem.data(), conn.colent2elem.size());
//       writearray2file(pde.datapath + "/cgent2dgent.bin", conn.cgent2dgent.data(), conn.cgent2dgent.size());
//     }        
    
    std::cout << "Finished buildConn.\n";
}


// int facepartition(int** elemcon, int** facecon, int** bf, int** f2t, int** facepartpts, int** facepartbnd, 
//                   const int* t, const int* f, const int* elempart, const int* localfaces,  
//                   const int* perm, const int* permind,  const int* bcm, int dim, int elemtype, 
//                   int nve, int nvf, int nfe, int npf, int npe, int ne, int ne1, int nbcm) 
// {
//       *bf = (int*)malloc(nfe * ne * sizeof(int));
//   
//       int* fi = (int*)malloc(nfe * ne * sizeof(int));      
//       select_columns(fi, f, elempart, nfe, ne);       
//       apply_bcm(*bf, fi, bcm, nve*ne, nbcm);       
//       CPUFREE(fi);
//                                           
//       int* ti = (int*)malloc(nve * ne * sizeof(int));                  
//       select_columns(ti, t, elempart, nve, ne); 
//       
//       int* f2t_tmp = (int*)malloc(4 * nfe * ne * sizeof(int));
//       int nf1 = mkf2e(f2t_tmp, ti, localfaces, ne, nve, nvf, nfe);       
//       int* ind = (int*)malloc(nf1 * sizeof(int));
//       int nf = 0;
//       for (int i=0; i<nf1; i++)
//         if (f2t_tmp[0 + 4*i] < ne1)
//           ind[nf++] = i;                
//       
//       *f2t = (int*)malloc(4 * nf * sizeof(int));
//       select_columns(*f2t, f2t_tmp, ind, 4, nf);        
//       CPUFREE(f2t_tmp);
//       
//       int nifaces = find(*f2t, 0, 4, nf, 2, 2);
//       int nbfaces = find(*f2t, -1, 4, nf, 2, 0);      
//       int* ifaces = (int*)malloc(nifaces * sizeof(int));
//       int* bfaces = (int*)malloc(nbfaces * sizeof(int));
//       find(ifaces, *f2t, 0, 4, nf, 2, 2);  // interior faces
//       find(bfaces, *f2t, -1, 4, nf, 2, 0); // boundary faces
//             
//       for (int i=0; i<nbfaces; i++) {
//         int e = *f2t[0 + 4*bfaces[i]];
//         int l = *f2t[1 + 4*bfaces[i]];
//         bfaces[i] = *bf[l + nfe*e];
//       }                    
//       
//       int* sortedbfaces = (int*)malloc(nbfaces * sizeof(int));
//       simple_bubble_sort(sortedbfaces, ind, bfaces, nbfaces); 
//       
//       *facepartbnd = (int*)malloc( (nbcm+1) * sizeof(int));
//       *facepartpts = (int*)malloc( (nbcm+1) * sizeof(int));
//       *facepartbnd[0] = 0;
//       *facepartpts[0] = nifaces;
//       int nbc = unique_count(&*facepartbnd[1], &*facepartpts[1], sortedbfaces, nbfaces);      
//             
//       int* bind = (int*)malloc(nf * sizeof(int));
//       for (int i=0; i<nifaces; i++) bind[i] = ifaces[i];
//       for (int i=0; i<nbfaces; i++) bind[nifaces+i] = bfaces[ind[i]];
//       permute_columns(*f2t, bind, 4, nf); 
//       
//       CPUFREE(ifaces); CPUFREE(bfaces); CPUFREE(ind); CPUFREE(bind);  CPUFREE(sortedbfaces);                
//       
//       fix_f2t_consistency(*f2t, elempart, nf); 
//       
//       *elemcon = (int*)malloc(npf * nfe * ne * sizeof(int)); 
//       *facecon = (int*)malloc(npf * 2 * nf * sizeof(int)); 
//             
//       build_connectivity(*elemcon, *facecon, *f2t, localfaces, ti, perm, 
//               permind, dim, elemtype, npf, nfe, npe, ne, nf); 
//                   
//       CPUFREE(ti);
//       
//       return nf;      
// }

#endif

