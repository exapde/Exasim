/**
 * @file makemesh.cpp
 * @brief Mesh generation and boundary/periodic face assignment utilities for finite element methods.
 *
 * This file contains functions for:
 *   - Determining local node indices for element faces (getelemface)
 *   - Computing face node coordinates and connectivity (compute_pf)
 *   - Evaluating expressions at mesh points (eval_expr)
 *   - Assigning boundary faces based on geometric expressions (assignboundaryfaces)
 *   - Sorting and comparing face node arrays (sort_facenodes, equal_faces)
 *   - Building face-to-element and element-to-element connectivity (mkf2e, mke2e)
 *   - Identifying and marking boundary faces (setboundaryfaces)
 *   - Handling periodic boundary conditions and matching periodic faces (setperiodicfaces)
 *   - Extracting interface elements for coupled interfaces (interface_elements)
 *   - Converting vector<string> to char** for boundary expressions (assignVectorToCharArray, freeCharArray)
 *   - Computing DG node coordinates from mesh and shape functions (compute_dgnodes)
 *   - Projecting DG nodes onto curved boundaries using distance functions (project_dgnodes_onto_curved_boundaries)
 *   - Initializing mesh structure from input parameters and PDE settings (initializeMesh)
 *
 * All arrays are assumed to be column-major and use zero-based indexing.
 * The code supports simplex and tensor elements in 1D, 2D, and 3D.
 * Boundary and periodic faces are assigned using user-provided expressions.
 * Curved boundaries are handled via projection using distance functions.
 *
 * Dependencies:
 *   - TinyExpr or similar expression parser (te_parser)
 *   - Mesh and PDE data structures
 *   - Utility functions for reading mesh and field data
 *
 * @author Cuong Nguyen
 * @date 2024
 */

#ifndef __MAKEMESH
#define __MAKEMESH

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
// void getelemface(int* face, int dim, int elemtype)
// {
//     int nvf = (dim == 3) ? (dim + elemtype) : dim;     
//     int nfe = dim + (dim-1)*elemtype + 1;    
// 
//     if (dim == 1) {
//         //face = (int*) malloc(sizeof(int) * (*nfe) * (nvf));
//         face[0] = 1;  // face 1 → node 1
//         face[1] = 2;  // face 2 → node 2
//     }
//     else if (dim == 2) {
//         if (elemtype == 0) {  // triangle
//             //face = (int*) malloc(sizeof(int) * (*nfe) * (nvf));
//             // Each column corresponds to a face; column-major layout
//             face[0 + 0*(nvf)] = 2;
//             face[1 + 0*(nvf)] = 3;
// 
//             face[0 + 1*(nvf)] = 3;
//             face[1 + 1*(nvf)] = 1;
// 
//             face[0 + 2*(nvf)] = 1;
//             face[1 + 2*(nvf)] = 2;
//         }
//         else if (elemtype == 1) {  // quadrilateral
//             //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
//             face[0 + 0*(nvf)] = 1;
//             face[1 + 0*(nvf)] = 2;
// 
//             face[0 + 1*(nvf)] = 2;
//             face[1 + 1*(nvf)] = 3;
// 
//             face[0 + 2*(nvf)] = 3;
//             face[1 + 2*(nvf)] = 4;
// 
//             face[0 + 3*(nvf)] = 4;
//             face[1 + 3*(nvf)] = 1;
//         }
//     }
//     else if (dim == 3) {
//         if (elemtype == 0) {  // tetrahedron
//             //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
//             face[0 + 0*(nvf)] = 2;
//             face[1 + 0*(nvf)] = 3;
//             face[2 + 0*(nvf)] = 4;
// 
//             face[0 + 1*(nvf)] = 1;
//             face[1 + 1*(nvf)] = 4;
//             face[2 + 1*(nvf)] = 3;
// 
//             face[0 + 2*(nvf)] = 1;
//             face[1 + 2*(nvf)] = 2;
//             face[2 + 2*(nvf)] = 4;
// 
//             face[0 + 3*(nvf)] = 1;
//             face[1 + 3*(nvf)] = 3;
//             face[2 + 3*(nvf)] = 2;
//         }
//         else if (elemtype == 1) {  // hexahedron
//             //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
//             face[0 + 0*(nvf)] = 1;
//             face[1 + 0*(nvf)] = 4;
//             face[2 + 0*(nvf)] = 3;
//             face[3 + 0*(nvf)] = 2;
// 
//             face[0 + 1*(nvf)] = 5;
//             face[1 + 1*(nvf)] = 6;
//             face[2 + 1*(nvf)] = 7;
//             face[3 + 1*(nvf)] = 8;
// 
//             face[0 + 2*(nvf)] = 1;
//             face[1 + 2*(nvf)] = 2;
//             face[2 + 2*(nvf)] = 6;
//             face[3 + 2*(nvf)] = 5;
// 
//             face[0 + 3*(nvf)] = 3;
//             face[1 + 3*(nvf)] = 4;
//             face[2 + 3*(nvf)] = 8;
//             face[3 + 3*(nvf)] = 7;
// 
//             face[0 + 4*(nvf)] = 2;
//             face[1 + 4*(nvf)] = 3;
//             face[2 + 4*(nvf)] = 7;
//             face[3 + 4*(nvf)] = 6;
// 
//             face[0 + 5*(nvf)] = 4;
//             face[1 + 5*(nvf)] = 1;
//             face[2 + 5*(nvf)] = 5;
//             face[3 + 5*(nvf)] = 8;
//         }
//     }
//     else {
//         fprintf(stderr, "Error: Only dim = 1, 2, or 3 supported.\n");
//         exit(1);
//     }
// 
//     for (int i=0; i<nvf*nfe; i++) face[i] -= 1;
// }

// Inputs:
// - int* t:       [nve * ne], element-to-node connectivity (column-major)
// - double* p:    [dim * np], coordinates of all points (column-major)
// - int* localfaces: [nvf * nfe], zero-based local face connectivity
// Outputs:
// - int* t2lf:    [nvf * nfe * ne], allocated by caller
// - double* pf:   [dim * nvf * nfe * ne], allocated by caller

void compute_pf(int* t2lf, double* pf,        
    const int* t, const double* p, const int* localfaces, 
    int nve, int nvf, int nfe, int dim, int ne   
) {
    // Step 1: Build t2lf: [nvf x nfe*ne]
    for (int e = 0; e < ne; ++e) {
        for (int f = 0; f < nfe; ++f) {
            for (int v = 0; v < nvf; ++v) {
                int lf = localfaces[v + nvf * f];      // zero-based local face node index
                int node = t[lf + nve * e];            // node index from element e
                int idx = v + nvf * (f + nfe * e);     // linear index in t2lf
                t2lf[idx] = node;
            }
        }
    }

    // Step 2: Build pf: [dim x nvf x nfe x ne]
    for (int e = 0; e < ne; ++e) {
        for (int f = 0; f < nfe; ++f) {
            for (int v = 0; v < nvf; ++v) {
                int node = t2lf[v + nvf * (f + nfe * e)];
                for (int d = 0; d < dim; ++d) {
                    int idx_pf = d + dim * (v + nvf * (f + nfe * e));
                    pf[idx_pf] = p[d + dim * node];
                }
            }
        }
    }
}

// Helper to evaluate expression string at a point x[dim]
bool eval_expr(const char* expr_str, const double* x, int dim) {
    double x0 = x[0], x1 = (dim > 1 ? x[1] : 0), x2 = (dim > 2 ? x[2] : 0);
    te_parser tep;
    tep.set_variables_and_functions({{"x", &x0}, {"y", &x1}, {"z", &x2}});
    
    /* This will compile the expression and check for errors. */
    auto result = tep.evaluate(expr_str);
    
    if (!tep.success()) {
        /* Show the user where the error is at. */
        std::cout<<expr_str<<std::endl;
        std::cout << "\t " << std::setfill(' ') <<
        std::setw(tep.get_last_error_position()) << '^' <<"\tError near here\n";      
        return EXIT_FAILURE;
    }
      
    return result;
    
//     te_variable vars[] = {
//         {"x", &x0}, {"y", &x1}, {"z", &x2}
//     };
//     int err;
//     te_expr* expr = te_compile(expr_str, vars, dim, &err);
//     if (!expr) {
//         fprintf(stderr, "TinyExpr error at %d in expression: %s\n", err, expr_str);
//         return false;
//     }
//     double result = te_eval(expr);
//     te_free(expr);
//     return result != 0.0;
}

void assignboundaryfaces(
    int* f,                     // [nfe x ne], zero-based, must be initialized to -1
    const int* f2t,             // [4 x nf], face-to-element connectivity
    const int* ind,             // [nb], indices of boundary faces   
    const double* pf,           // [dim x nvf x nfe x ne], column-major
    char** bndexpr,       // boundary expressions, e.g., {"x*x + y*y <= 1"}
    int dim, int nvf, int nfe, int ne, int nb, int nbndexpr) 
{
    for (int i = 0; i < nfe * ne; ++i)
        f[i] = -1;

    int ind2 = 1;
    if (nvf==1) ind2 = 0;
            
    for (int ii = 0; ii < nb; ++ii) {
        int b = ind[ii];
        int e = f2t[0 + 4 * b];
        int l = f2t[1 + 4 * b];

        for (int k = 0; k < nbndexpr; ++k) {
            bool ok = true;
            int check_pts[3] = {0, ind2, nvf - 1};

            for (int m = 0; m < 3; ++m) {
                int i = check_pts[m];
                // pf: dim x nvf x nfe x ne
                const double* x = &pf[dim * (i + nvf * (l + nfe * e))];
                //printf("%d %d %d %d %g %g %s\n", nvf, e, l, i, x[0], x[1], bndexpr[k]);
                //printf("%d %d %d %d %g %s\n", nvf, e, l, i, x[0], bndexpr[k]);
                if (!eval_expr(bndexpr[k], x, dim)) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                f[l + nfe * e] = k;
                break;
            }
        }
    }
}

void assignboundaryfaces(
    int* f,                     // [nfe x ne], zero-based, must be initialized to -1
    const int* t,               // [nve x ne], element-to-node connectivity
    const int* f2t,             // [4 x nf], face-to-element connectivity
    const int* ind,             // [nb], indices of boundary faces   
    const int* localfaces,      // [nvf x nfe], local face connectivity
    const double* p,           // [dim x np]
    char** bndexpr,       // boundary expressions, e.g., {"x*x + y*y <= 1"}
    int dim, int nve, int nvf, int nfe, int ne, int nb, int nbndexpr) 
{
    for (int i = 0; i < nfe * ne; ++i)
        f[i] = -1;

    int ind2 = 1;
    if (nvf==1) ind2 = 0;
            
    for (int ii = 0; ii < nb; ++ii) {
        int b = ind[ii];
        int e = f2t[0 + 4 * b];
        int l = f2t[1 + 4 * b];

        for (int k = 0; k < nbndexpr; ++k) {
            bool ok = true;
            int check_pts[3] = {0, ind2, nvf - 1};

            for (int m = 0; m < 3; ++m) {
                int i = check_pts[m];

                // Get the local node index for this face
                int local_node = localfaces[i + l * nvf];
                int node = t[local_node + nve * e];   // node index from element e
                const double* x = &p[dim * node];

                if (!eval_expr(bndexpr[k], x, dim)) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                f[l + nfe * e] = k;
                break;
            }
        }
    }
}

// bool equal_faces(const int* a, const int* b, int nnf) {
//     for (int i = 0; i < nnf; i++) {
//         if (a[i] != b[i]) return false;
//     }
//     return true;
// }
// 
// // Sort small arrays (nnf elements) using bubble sort
// void sort_facenodes(int* face, int nnf) {
//     for (int i = 0; i < nnf - 1; ++i) {
//         for (int j = 0; j < nnf - i - 1; ++j) {
//             if (face[j] > face[j + 1]) {
//                 int tmp = face[j];
//                 face[j] = face[j + 1];
//                 face[j + 1] = tmp;
//             }
//         }
//     }
// }

// int mkf2e(int* f2e, const int* e2n, const int* local_faces, int ne, int nne, int nnf, int nfe) 
// {
//     int max_faces = ne * nfe;
//     int* face_nodes = (int*)malloc(sizeof(int) * nnf * max_faces); // stores sorted face nodes
//     int count = 0;
// 
//     for (int elem = 0; elem < ne; ++elem) {
//         for (int lf = 0; lf < nfe; ++lf) {
//             int temp_face[8];  
//             for (int i = 0; i < nnf; ++i) {
//                 int local_node = local_faces[i + lf * nnf];
//                 temp_face[i] = e2n[elem * nne + local_node];
//             }
// 
//             sort_facenodes(temp_face, nnf);
// 
//             // Check if face already exists
//             int match = -1;
//             for (int f = 0; f < count; ++f) {
//                 if (equal_faces(temp_face, &face_nodes[f * nnf], nnf)) {
//                     match = f;
//                     break;
//                 }
//             }
// 
//             if (match >= 0) {
//                 // Second element sharing the face
//                 f2e[2 + match * 4] = elem;  // element on other side
//                 f2e[3 + match * 4] = lf;    // local face index in that element
//             } else {
//                 // New face
//                 for (int i = 0; i < nnf; ++i)
//                     face_nodes[count * nnf + i] = temp_face[i];
// 
//                 f2e[0 + count * 4] = elem;  // first element
//                 f2e[1 + count * 4] = lf;    // local face index in first element
//                 f2e[2 + count * 4] = -1;    // no neighbor yet
//                 f2e[3 + count * 4] = -1;
//                 ++count;
//             }
//         }
//     }
// 
//     CPUFREE(face_nodes);
//     return count;    
// }
// 
// int mke2e(int* e2e, const int* e2n, const int* local_faces, int ne, int nne, int nnf, int nfe) 
// {
//     int max_faces = ne * nfe;
//     int* face_nodes = (int*)malloc(sizeof(int) * nnf * max_faces);
//     int* face_owner = (int*)malloc(sizeof(int) * max_faces);
//     int* face_lf = (int*)malloc(sizeof(int) * max_faces);
//     int face_count = 0;
// 
//     for (int elem = 0; elem < ne; ++elem) {
//         for (int lf = 0; lf < nfe; ++lf) {
//             int temp_face[8];
//             for (int i = 0; i < nnf; ++i) {
//                 int local_node = local_faces[i + lf * nnf];
//                 temp_face[i] = e2n[elem * nne + local_node];
//             }
// 
//             sort_facenodes(temp_face, nnf);
// 
//             int match = -1;
//             for (int f = 0; f < face_count; ++f) {
//                 if (equal_faces(temp_face, &face_nodes[f * nnf], nnf)) {
//                     match = f;
//                     break;
//                 }
//             }
// 
//             if (match >= 0) {
//                 // Face shared with another element
//                 int e1 = face_owner[match];
//                 int l1 = face_lf[match];
//                 int e2 = elem;
//                 int l2 = lf;
// 
//                 e2e[e1 * nfe + l1] = e2;
//                 e2e[e2 * nfe + l2] = e1;
//             } else {
//                 // New face
//                 for (int i = 0; i < nnf; ++i)
//                     face_nodes[face_count * nnf + i] = temp_face[i];
//                 face_owner[face_count] = elem;
//                 face_lf[face_count] = lf;
//                 e2e[elem * nfe + lf] = -1;  // boundary face
//                 face_count++;
//             }
//         }
//     }
// 
//     CPUFREE(face_nodes);
//     CPUFREE(face_owner);
//     CPUFREE(face_lf);
//     return face_count;
// }

// Hashable key of fixed max size 8 (adjust if needed)
// struct FaceKey {
//     std::array<int,8> a{};
//     uint8_t len{}; // nnf
// 
//     bool operator==(const FaceKey& o) const noexcept {
//         if (len != o.len) return false;
//         for (uint8_t i = 0; i < len; ++i) if (a[i] != o.a[i]) return false;
//         return true;
//     }
// };
// struct FaceKeyHash {
//     size_t operator()(const FaceKey& k) const noexcept {
//         // simple 64-bit mix across len entries
//         uint64_t h = 0x9E3779B97F4A7C15ull ^ k.len;
//         for (uint8_t i = 0; i < k.len; ++i) {
//             uint64_t x = static_cast<uint64_t>(k.a[i]) + 0x9E3779B97F4A7C15ull;
//             x ^= (x >> 30); x *= 0xBF58476D1CE4E5B9ull;
//             x ^= (x >> 27); x *= 0x94D049BB133111EBull;
//             x ^= (x >> 31);
//             h ^= x + 0x9E3779B97F4A7C15ull + (h<<6) + (h>>2);
//         }
//         return static_cast<size_t>(h);
//     }
// };
// 
// int mkf2e_hash(int* f2e,
//                const int* e2n, const int* local_faces,
//                int ne, int nne, int nnf, int nfe)
// {
//     using Val = std::pair<int,int>; // (elem, lf)
//     std::unordered_map<FaceKey, Val, FaceKeyHash> map;
//     map.reserve(static_cast<size_t>(ne) * nfe * 1.3);
// 
//     int out = 0;
// 
//     auto make_key = [&](int e, int lf) {
//         FaceKey k; k.len = static_cast<uint8_t>(nnf);
//         // gather face nodes
//         for (int i = 0; i < nnf; ++i) {
//             int ln = local_faces[lf*nnf + i];
//             k.a[i] = e2n[e*nne + ln];
//         }
//         std::sort(k.a.begin(), k.a.begin() + nnf);
//         return k;
//     };
// 
//     for (int e = 0; e < ne; ++e) {
//         for (int lf = 0; lf < nfe; ++lf) {
//             FaceKey k = make_key(e, lf);
//             auto it = map.find(k);
//             if (it == map.end()) {
//                 map.emplace(std::move(k), Val{e, lf});
//             } else {
//                 // found the interior faces
//                 const auto [e0, lf0] = it->second;
//                 f2e[4*out + 0] = e0;  f2e[4*out + 1] = lf0;
//                 f2e[4*out + 2] = e;   f2e[4*out + 3] = lf;
//                 map.erase(it);
//                 ++out;
//             }
//         }
//     }
// 
//     // Remaining entries in map are boundary faces
//     for (const auto& kv : map) {
//         const auto [e0, lf0] = kv.second;
//         f2e[4*out + 0] = e0;  f2e[4*out + 1] = lf0;
//         f2e[4*out + 2] = -1;  f2e[4*out + 3] = -1;
//         ++out;
//     }
// 
//     return out;
// }
// 
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

int setboundaryfaces(int* f, int* t2lf, int *face, const double* p, const int* t,                    
    char** bndexpr, int dim, int elemtype, int ne, int nbndexpr) 
{
    int nve = (elemtype==0) ? (dim + 1) : std::pow(2, dim);    
    int nvf = (dim == 3) ? (dim + elemtype) : dim;     
    int nfe = dim + (dim-1)*elemtype + 1;
    
    //int* face = (int*) malloc(sizeof(int) * nfe * nvf);
    getelemface(face, dim, elemtype);    
    
    for (int e = 0; e < ne; ++e) {
        for (int f = 0; f < nfe; ++f) {
            for (int v = 0; v < nvf; ++v) {
                int lf = face[v + nvf * f];      // zero-based local face node index
                int node = t[lf + nve * e];            // node index from element e
                int idx = v + nvf * (f + nfe * e);     // linear index in t2lf
                t2lf[idx] = node;
            }
        }
    }    

    int* f2e = (int*) malloc(sizeof(int) * 4 * nfe * ne);
    int nf = mkf2e_hash(f2e, t, face, ne, nve, nvf, nfe); 
    
    int* ind = (int*) malloc(sizeof(int) * nf);
    int nb = 0;
    for (int i=0; i<nf; i++) ind[i] = -1;
    for (int i=0; i<nf; i++) 
      if (f2e[2 + 4*i] == -1) ind[nb++] = i;

    assignboundaryfaces(f, t, f2e, ind, face, p, bndexpr, dim, nve, nvf, nfe, ne, nb, nbndexpr);    

    CPUFREE(f2e);
    CPUFREE(ind);
    
    return nf;
}

void setperiodicfaces(int* f, int* t, const double* p, const int* t2fl, 
    const int* prd_f1, const int* prd_f2, char** expr1, char** expr2,    
    int dim, int elemtype, int np, int ne, int nprd, int ncomp) 
{  
    int nve = (elemtype==0) ? (dim + 1) : std::pow(2, dim);    
    int nvf = (dim == 3) ? (dim + elemtype) : dim;     
    int nfe = dim + (dim-1)*elemtype + 1;
  
    int nface = nfe * ne;

    for (int i = 0; i < nprd; ++i) {
        // Collect indices of faces on the first and second periodic boundary
        int* i1 = (int*)malloc(nface * sizeof(int));
        int* i2 = (int*)malloc(nface * sizeof(int));
        int n1 = 0, n2 = 0;
        for (int j = 0; j < nface; ++j) {
            if (f[j] == prd_f1[i]-1) {
                f[j] = -f[j]-2;
                i1[n1++] = j;
            } else if (f[j] == prd_f2[i]-1) {
                f[j] = -f[j]-2;
                i2[n2++] = j;
            }
        }

        // Get vertices v1 and v2 from t2fl
        int* v1 = (int*)malloc(n1 * nvf * sizeof(int));
        int* v2 = (int*)malloc(n2 * nvf * sizeof(int));
        int v1count = 0, v2count = 0;
        for (int j = 0; j < n1; ++j)
            for (int k = 0; k < nvf; ++k)
                v1[v1count++] = t2fl[k + nvf * i1[j]];
        for (int j = 0; j < n2; ++j)
            for (int k = 0; k < nvf; ++k)
                v2[v2count++] = t2fl[k + nvf * i2[j]];

        int v1len = unique_ints(v1, v1count);
        int v2len = unique_ints(v2, v2count);

        // Evaluate expressions on p(v1) and p(v2)
        double* p1 = (double*)malloc(v1len * dim * sizeof(double));
        double* p2 = (double*)malloc(v2len * dim * sizeof(double));
        for (int j = 0; j < v1len; ++j)
            for (int d = 0; d < dim; ++d)
                p1[j * dim + d] = p[d + dim * v1[j]];
        for (int j = 0; j < v2len; ++j)
            for (int d = 0; d < dim; ++d)
                p2[j * dim + d] = p[d + dim * v2[j]];

        double x=0, y=0, z=0;
        te_parser tep;        
        tep.set_variables_and_functions({{"x", &x}, {"y", &y}, {"z", &z}});
                
        double* q1 = (double*)malloc(v1len * ncomp * sizeof(double));
        double* q2 = (double*)malloc(v2len * ncomp * sizeof(double));

        for (int c = 0; c < ncomp; ++c) {
            auto result = tep.evaluate(expr1[i * ncomp + c]);
            if (!tep.success()) {
                /* Show the user where the error is at. */
                std::cout << "\t " << std::setfill(' ') <<
                std::setw(tep.get_last_error_position()) << '^' <<"\tError near here\n";      
                error("TinyExpr Failure");
            }            
            for (int j = 0; j < v1len; ++j) {
                x = (dim > 0) ? p1[j * dim + 0] : 0;
                y = (dim > 1) ? p1[j * dim + 1] : 0;
                z = (dim > 2) ? p1[j * dim + 2] : 0;
                q1[j * ncomp + c] = tep.evaluate(expr1[i * ncomp + c]);            
                
            }
            result = tep.evaluate(expr2[i * ncomp + c]);
            if (!tep.success()) {
                /* Show the user where the error is at. */
                std::cout << "\t " << std::setfill(' ') <<
                std::setw(tep.get_last_error_position()) << '^' <<"\tError near here\n";      
                error("TinyExpr Failure");
            }                        
            for (int j = 0; j < v2len; ++j) {
                x = (dim > 0) ? p2[j * dim + 0] : 0;
                y = (dim > 1) ? p2[j * dim + 1] : 0;
                z = (dim > 2) ? p2[j * dim + 2] : 0;
                q2[j * ncomp + c] = tep.evaluate(expr2[i * ncomp + c]);                      
            }
        }
        
//     te_parser tep;
//     tep.set_variables_and_functions({{"x", &x0}, {"y", &x1}, {"z", &x2});
//     
//     /* This will compile the expression and check for errors. */
//     auto result = tep.evaluate(expr_str);
//     
//     if (!tep.success()) {
//         /* Show the user where the error is at. */
//         std::cout << "\t " << std::setfill(' ') <<
//         std::setw(tep.get_last_error_position()) << '^' <<"\tError near here\n";      
//         return EXIT_FAILURE;
//     }
//       
//     return result;
        
//         double x, y, z;
//         te_variable vars[] = { {"x", &x}, {"y", &y}, {"z", &z} };
//         int len = (dim >= 3) ? 3 : dim;
// 
//         te_expr* ex1 = te_compile(expr1[i], vars, len, NULL);
//         te_expr* ex2 = te_compile(expr2[i], vars, len, NULL);
// 
//         if (!ex1 || !ex2) {
//             fprintf(stderr, "TinyExpr parse error in periodic expression %d\n", i);
//             exit(1);
//         }        
//         
//         double* q1 = (double*)malloc(v1len * ncomp * sizeof(double));
//         double* q2 = (double*)malloc(v2len * ncomp * sizeof(double));
// 
//         for (int c = 0; c < ncomp; ++c) {
//             te_expr* ex1 = te_compile(expr1[i * ncomp + c], vars, dim, NULL);
//             te_expr* ex2 = te_compile(expr2[i * ncomp + c], vars, dim, NULL);
//             if (!ex1 || !ex2) {
//                 fprintf(stderr, "TinyExpr parse error in periodic expression %d component %d\n", i, c);
//                 exit(1);
//             }
//             for (int j = 0; j < v1len; ++j) {
//                 x = (dim > 0) ? p1[j * dim + 0] : 0;
//                 y = (dim > 1) ? p1[j * dim + 1] : 0;
//                 z = (dim > 2) ? p1[j * dim + 2] : 0;
//                 q1[j * ncomp + c] = te_eval(ex1);
//             }
//             for (int j = 0; j < v2len; ++j) {
//                 x = (dim > 0) ? p2[j * dim + 0] : 0;
//                 y = (dim > 1) ? p2[j * dim + 1] : 0;
//                 z = (dim > 2) ? p2[j * dim + 2] : 0;
//                 q2[j * ncomp + c] = te_eval(ex2);
//             }
//             te_free(ex1);
//             te_free(ex2);
//         }
        
//         double* q1 = (double*)malloc(v1len * sizeof(double));
//         double* q2 = (double*)malloc(v2len * sizeof(double));
//         
//         for (int j = 0; j < v1len; ++j) {
//             x = (dim > 0) ? p1[j * dim + 0] : 0;
//             y = (dim > 1) ? p1[j * dim + 1] : 0;
//             z = (dim > 2) ? p1[j * dim + 2] : 0;
//             q1[j] = te_eval(ex1); // may update values
//         }
//         for (int j = 0; j < v2len; ++j) {
//             x = (dim > 0) ? p2[j * dim + 0] : 0;
//             y = (dim > 1) ? p2[j * dim + 1] : 0;
//             z = (dim > 2) ? p2[j * dim + 2] : 0;
//             q2[j] = te_eval(ex2);
//         }

        // Match and permute vertices
        int* in = (int*)malloc(v1len * sizeof(int));
        xiny<double>(in, q1, q2, v1len, v2len, ncomp);
        for (int j = 0; j < v1len; ++j) {
          if (in[j]>=0) {
            int old_val = v2[in[j]];
            int new_val = v1[j];
            for (int k = 0; k < nve * ne; ++k) {
                if (t[k] == old_val) t[k] = new_val;
            }
          }
        }

//         te_free(ex1);
//         te_free(ex2);

        CPUFREE(i1); CPUFREE(i2); CPUFREE(v1); CPUFREE(v2); CPUFREE(p1); CPUFREE(p2); CPUFREE(q1); CPUFREE(q2); CPUFREE(in);
    }
}

void interface_elements(vector<int>& inte, vector<int>& intl, vector<int>& f, int nfe, int ne, int coupledinterface)
{
    int n=0;
    for (int i=0; i<ne; i++)
      for (int j=0; j<nfe; j++)
        if (f[j + nfe*i] == coupledinterface) n += 1;

    inte.resize(n);
    intl.resize(n);
    n = 0; 
    for (int i=0; i<ne; i++)
      for (int j=0; j<nfe; j++) 
        if (f[j + nfe*i] == coupledinterface) {
          inte[n] = i;
          intl[n] = j;
          n += 1;
        }    
}

// Function to convert vector<string> A to char** B
void assignVectorToCharArray(const std::vector<std::string>& A, char*** B) {
    *B = new char*[A.size()];
    for (size_t i = 0; i < A.size(); ++i) {
        (*B)[i] = new char[A[i].size() + 1];  // +1 for null terminator
        std::strcpy((*B)[i], A[i].c_str());
    }
}

// Optional: helper function to free the allocated memory
void freeCharArray(char** B, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        delete[] B[i];
    }
    delete[] B;
}

// All arrays are assumed to be column-major, 0-based indexing
// p    : double[nd * npv]            -- flattened 2D array, p[dim + nd * j]
// t    : int[nve * ne]               -- connectivity: t[node + nve * e]
// philocal: double[npl * nve]        -- shape function values: phi[i + npl * node]
// dgnodes : double[npl * nd * ne]    -- output array
void compute_dgnodes(double* dgnodes, const double* p, const int* t, 
                     const double* philocal, int npl, int nd, int ne, int nve)
{
    // Zero out dgnodes
    for (int i = 0; i < npl * nd * ne; ++i)
        dgnodes[i] = 0.0;
    
    for (int e = 0; e < ne; ++e) {
        for (int dim = 0; dim < nd; ++dim) {
            for (int i = 0; i < npl; ++i) {
              for (int node = 0; node < nve; ++node) {             
                  int global_node = t[node + nve * e];                     
                  int idx = i + npl * (dim + nd * e);   // dgnodes(i,dim,e)
                  dgnodes[idx] += philocal[i + npl * node] * p[dim + nd * global_node];
              }
            }
        }
    }        
}

void project_dgnodes_onto_curved_boundaries(double* dgnodes, const int* f, const int* perm, const int* curvedboundary,
                         char** fd_exprs, int nd, int porder, int npe, int npf, int nfe, int ne) 
{
    if (porder <= 1) return;

    int has_curved = 0;
    for (int k = 0; k < nfe * ne; ++k)
        if (f[k] > -1 && curvedboundary[f[k]] != 0)
            has_curved = 1;

    if (!has_curved) return;

    // printf("Project dgnodes onto the curved boundaries...\n");

    double eps = 1e-14;
    double deps;

    for (int e = 0; e < ne; ++e) {
        for (int j = 0; j < nfe; ++j) {
            int fid = f[j + nfe * e];
            if (fid < 0) continue;

            int k = fid;
            if (curvedboundary[k] != 1) continue;

            const char* expr_str = fd_exprs[k];
            double x=0, y=0, z=0;
            
            te_parser tep;
            tep.set_variables_and_functions({{"x", &x}, {"y", &y}, {"z", &z}});
            
            /* This will compile the expression and check for errors. */
            auto result = tep.evaluate(expr_str);
            
            if (!tep.success()) {
                /* Show the user where the error is at. */
                std::cout << "\t " << std::setfill(' ') <<
                std::setw(tep.get_last_error_position()) << '^' <<"\tError near here\n";      
                error("TinyExpr Failure");
            }
            
//             te_variable vars[] = { {"x", &x}, {"y", &y}, {"z", &z} };
//             int err;
//             te_expr* expr = te_compile(expr_str, vars, nd, &err);
//             if (!expr) {
//                 fprintf(stderr, "TinyExpr compile error at %d for expression: %s\n", err, expr_str);
//                 continue;
//             }

            // collect coordinates: p[npts][nd]
            int npts = npf;
            double* p = (double*) malloc(sizeof(double) * npts * nd);
            for (int i = 0; i < npts; ++i) {
                int idx = perm[i + npf * j];
                for (int d = 0; d < nd; ++d)
                    p[i + npts * d] = dgnodes[idx + npe * (d + nd * e)];
            }

            // estimate deps
            double pmax = -1e20, pmin = 1e20;
            for (int i = 0; i < npts * nd; ++i) {
                if (p[i] > pmax) pmax = p[i];
                if (p[i] < pmin) pmin = p[i];
            }
            deps = sqrt(eps) * (pmax - pmin);

            // evaluate distance function 
            double* d = (double*) malloc(sizeof(double) * npts);
            for (int i = 0; i < npts; ++i) {
                x = (nd > 0) ? p[i + npts * 0] : 0;
                y = (nd > 1) ? p[i + npts * 1] : 0;
                z = (nd > 2) ? p[i + npts * 2] : 0;
                //d[i] = te_eval(expr);
                d[i] = tep.evaluate(expr_str);
            }

            // evaluate finite-difference gradients
            double* dgrad = (double*) malloc(sizeof(double) * npts * nd);
            for (int ddir = 0; ddir < nd; ++ddir) {
                for (int i = 0; i < npts; ++i)
                    p[i + npts * ddir] += deps;

                for (int i = 0; i < npts; ++i) {
                    x = (nd > 0) ? p[i + npts * 0] : 0;
                    y = (nd > 1) ? p[i + npts * 1] : 0;
                    z = (nd > 2) ? p[i + npts * 2] : 0;
                    //double dshift = te_eval(expr);
                    double dshift = tep.evaluate(expr_str);
                    dgrad[i + npts * ddir] = (dshift - d[i]) / deps;
                }

                for (int i = 0; i < npts; ++i)
                    p[i + npts * ddir] -= deps;
            }

            // project points onto the boundary
            for (int i = 0; i < npts; ++i) {
                double sumsq = 0;
                for (int ddir = 0; ddir < nd; ++ddir)
                    sumsq += dgrad[i + npts * ddir] * dgrad[i + npts * ddir];
                if (sumsq == 0) sumsq = 1;
                for (int ddir = 0; ddir < nd; ++ddir)
                    p[i + npts * ddir] -= (d[i] * dgrad[i + npts * ddir]) / sumsq;
            }

            // update dgnodes
            for (int i = 0; i < npts; ++i) {
                int idx = perm[i + npf * j];
                for (int ddir = 0; ddir < nd; ++ddir)
                    dgnodes[idx + npe * (ddir + nd * e)] = p[i + npts * ddir];
            }

            CPUFREE(d);
            CPUFREE(dgrad);
            CPUFREE(p);
        }
    }
    
    //std::cout << "Finished project_dgnodes_onto_curved_boundaries.\n";
}

Mesh initializeMesh(InputParams& params, PDE& pde)
{    
    Mesh mesh;
    readMeshFromFile(make_path(pde.datapath, pde.meshfile), mesh);                 
    if (pde.xdgfile != "") {
        std::cout << "Reading xdg file.\n";
        readFieldFromBinaryFile(make_path(pde.datapath, pde.xdgfile), mesh.xdg, mesh.xdgdims);
    }
    if (pde.udgfile != "") {
        std::cout << "Reading udg file.\n";
        readFieldFromBinaryFile(make_path(pde.datapath, pde.udgfile), mesh.udg, mesh.udgdims);
    }
    if (pde.vdgfile != "") {
        std::cout << "Reading vdg file.\n";
        readFieldFromBinaryFile(make_path(pde.datapath, pde.vdgfile), mesh.vdg, mesh.vdgdims);
    }
    if (pde.wdgfile != "") {
        std::cout << "Reading wdg file.\n";
        readFieldFromBinaryFile(make_path(pde.datapath, pde.wdgfile), mesh.wdg, mesh.wdgdims);    
    }
    if (pde.uhatfile != "") {
        std::cout << "Reading uhat file.\n";
        readFieldFromBinaryFile(make_path(pde.datapath, pde.uhatfile), mesh.uhat, mesh.uhatdims);    
    }
    if (pde.partitionfile != "") {
        std::cout << "Reading partition file.\n";
        readPartitionFromFile(make_path(pde.datapath, pde.partitionfile), mesh.elem2cpu, mesh.ne);    
    }
    
    mesh.boundaryConditions = params.boundaryConditions;
    mesh.curvedBoundaries = params.curvedBoundaries;
    mesh.periodicBoundaries1 = params.periodicBoundaries1;
    mesh.periodicBoundaries2 = params.periodicBoundaries2;
    mesh.cartGridPart = params.cartGridPart;
    mesh.interfaceConditions = params.interfaceConditions;
        
    assignVectorToCharArray(params.boundaryExprs, &mesh.boundaryExprs);
    assignVectorToCharArray(params.curvedBoundaryExprs, &mesh.curvedBoundaryExprs);
    assignVectorToCharArray(params.periodicExprs1, &mesh.periodicExprs1);
    assignVectorToCharArray(params.periodicExprs2, &mesh.periodicExprs2);
        
    mesh.nbndexpr = params.boundaryExprs.size();
    mesh.nbcm = params.boundaryConditions.size();
    mesh.nprdexpr = params.periodicBoundaries1.size();    
    mesh.nprdcom = (mesh.nprdexpr == 0) ? 0 : params.periodicExprs1.size()/mesh.nprdexpr;
    if (mesh.nbndexpr != mesh.nbcm) 
        error("boundaryconditions and boundaryexpressions are not the same size. Exiting.\n");
                    
    mesh.dim = mesh.nd;
    pde.nve = mesh.nve; pde.np = mesh.np; pde.ne = mesh.ne; pde.elemtype = mesh.elemtype;                
    pde.nd = mesh.dim; pde.ncx = mesh.dim;
    if (pde.model=="ModelC" || pde.model=="modelC") {
        pde.wave = 0;
        pde.nc = pde.ncu;
    } else if (pde.model=="ModelD" || pde.model=="modelD") {     
        pde.wave = 0;
        pde.nc = (pde.ncu)*(pde.nd+1);
    } else if (pde.model=="ModelW" || pde.model=="modelW") {
        pde.tdep = 1;
        pde.wave = 1;
        pde.nc = (pde.ncu)*(pde.nd+1);
    }
    pde.ncq = pde.nc - pde.ncu;
    pde.nch  = pde.ncu;               
    
    std::cout << "Finished reading mesh file.\n";

//     cout << "Read mesh with:" << endl;
//     cout << "  nd  = " << mesh.nd << endl;
//     cout << "  np  = " << mesh.np << endl;    
//     cout << "  nve = " << mesh.nve << endl;
//     cout << "  ne  = " << mesh.ne << endl;    
//     cout << "  nfe = " << mesh.nfe << endl;
//     cout << "  nvf = " << mesh.nvf << endl;
//     cout << "  elemtype = " << mesh.elemtype << endl;
//     
//     for (size_t i = 0; i < params.boundaryExprs.size(); ++i) {
//         std::cout << mesh.boundaryExprs[i] << std::endl;
//     }    
//     print2darray(mesh.p.data(), mesh.dim, mesh.np);
//     print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
    
//     mesh.f.resize(mesh.nfe * mesh.ne);
//     mesh.t2lf.resize(mesh.nvf * mesh.nfe * mesh.ne);
//     mesh.localfaces.resize(mesh.nvf * mesh.nfe);
//     
//     mesh.nf = setboundaryfaces(mesh.f.data(), mesh.t2lf.data(), mesh.localfaces.data(), mesh.p.data(), 
//                mesh.t.data(), mesh.boundaryExprs, mesh.dim, mesh.elemtype, mesh.ne, mesh.nbndexpr); 
//     
//     std::cout << "Finished setboundaryfaces.\n";
// 
//     //print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
//     
//     if (mesh.nprdexpr > 0) {
//         setperiodicfaces(mesh.f.data(), mesh.t.data(), mesh.p.data(), mesh.t2lf.data(), 
//             mesh.periodicBoundaries1.data(), mesh.periodicBoundaries2.data(), mesh.periodicExprs1, 
//             mesh.periodicExprs2, mesh.dim, mesh.elemtype, mesh.np, mesh.ne, mesh.nprdexpr, mesh.nprdcom); 
//         std::cout << "Finished setperiodicfaces.\n";
//     }
//     
//     //print2iarray(mesh.f.data(), mesh.nfe, mesh.ne);
//     //print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
//     
//     if (pde.coupledinterface>0) {
//         interface_elements(mesh.inte, mesh.intl, mesh.f, mesh.nfe, mesh.ne, pde.coupledinterface);
//         std::cout << "Finished interface_elements.\n";
//     }
//     
//     std::cout << "Finished initializeMesh.\n";
//     
    
    return mesh;
}


#endif

