#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cmath>

#include "tinyexpr.h"

template<typename T>
void xiny(int* out, const T* A, const T* B, int m, int n, int dim, dstype tol = 1e-12) {
    for (int i = 0; i < n; ++i) {
        out[i] = -1;
        for (int j = 0; j < m; ++j) {
            bool match = true;
            for (int d = 0; d < dim; ++d) {
                if (std::abs(B[i * dim + d] - A[j * dim + d]) > tol) {
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

int divide_interval(int* out_intervals, int n, int m) {
    if (n <= 0 || m <= 0) return 0;

    int num_intervals = (n + m - 1) / m; // ceil(n/m)

    int base = n / num_intervals;
    int rem = n % num_intervals;
    int current = 1;

    for (int i = 0; i < num_intervals; ++i) {
        int len = base + (i < rem ? 1 : 0);
        out_intervals[2 * i] = current;
        out_intervals[2 * i + 1] = current + len - 1;
        current += len;
    }

    return num_intervals;
}

int mkfaceblocks(int** out_nm, const int* mf, const int* bcm, int nmf_len, int ns) {
    if (ns <= 0) ns = 2048;  // default value

    int max_blocks = 0;
    for (int i = 0; i < nmf_len - 1; ++i)
        max_blocks += (mf[i + 1] - mf[i] + ns - 1) / ns;

    int* nm = (int*)malloc(3 * max_blocks * sizeof(int));
    int count = 0;

    for (int i = 0; i < nmf_len - 1; ++i) {
        int nf = mf[i + 1] - mf[i];
        int num_intervals = (nf + ns - 1) / ns; 
        int* intervals = (int*)malloc(2 * num_intervals * sizeof(int)); 
        int nblocks = divide_interval(intervals, nf, ns);

        for (int j = 0; j < nblocks; ++j) {
            int start = mf[i] + intervals[2 * j];       // 0-based
            int end   = mf[i] + intervals[2 * j + 1];    // 0-based
            nm[3 * count + 0] = start;
            nm[3 * count + 1] = end;
            nm[3 * count + 2] = bcm[i];  // boundary code
            count++;
        }

        free(intervals);
    }

    *out_nm = (int*)realloc(nm, 3 * count * sizeof(int));

    return count; 
}

void cumsum_int(const int* in, int* out, int n) {
    if (n <= 0) return;
    out[0] = in[0];
    for (int i = 1; i < n; ++i) {
        out[i] = out[i - 1] + in[i];
    }
}

// Sort and remove duplicates
int unique_ints(int* arr, int n) {
    // Simple insertion sort
    for (int i = 1; i < n; ++i) {
        int key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
    // Remove duplicates
    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (m == 0 || arr[i] != arr[m - 1]) {
            arr[m++] = arr[i];
        }
    }
    return m;
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
    //int nfe = dim + (dim-1)*elemtype + 1;
          
    if (dim == 1) {
        //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
        face[0] = 1;  // face 1 → node 1
        face[1] = 2;  // face 2 → node 2
    }
    else if (dim == 2) {
        if (elemtype == 0) {  // triangle
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            // Each column corresponds to a face; column-major layout
            face[0 + 0*(*nvf)] = 2;
            face[1 + 0*(*nvf)] = 3;

            face[0 + 1*(*nvf)] = 3;
            face[1 + 1*(*nvf)] = 1;

            face[0 + 2*(*nvf)] = 1;
            face[1 + 2*(*nvf)] = 2;
        }
        else if (elemtype == 1) {  // quadrilateral
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(*nvf)] = 1;
            face[1 + 0*(*nvf)] = 2;

            face[0 + 1*(*nvf)] = 2;
            face[1 + 1*(*nvf)] = 3;

            face[0 + 2*(*nvf)] = 3;
            face[1 + 2*(*nvf)] = 4;

            face[0 + 3*(*nvf)] = 4;
            face[1 + 3*(*nvf)] = 1;
        }
    }
    else if (dim == 3) {
        if (elemtype == 0) {  // tetrahedron
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(*nvf)] = 2;
            face[1 + 0*(*nvf)] = 3;
            face[2 + 0*(*nvf)] = 4;

            face[0 + 1*(*nvf)] = 1;
            face[1 + 1*(*nvf)] = 4;
            face[2 + 1*(*nvf)] = 3;

            face[0 + 2*(*nvf)] = 1;
            face[1 + 2*(*nvf)] = 2;
            face[2 + 2*(*nvf)] = 4;

            face[0 + 3*(*nvf)] = 1;
            face[1 + 3*(*nvf)] = 3;
            face[2 + 3*(*nvf)] = 2;
        }
        else if (elemtype == 1) {  // hexahedron
            //face = (int*) malloc(sizeof(int) * (*nfe) * (*nvf));
            face[0 + 0*(*nvf)] = 1;
            face[1 + 0*(*nvf)] = 4;
            face[2 + 0*(*nvf)] = 3;
            face[3 + 0*(*nvf)] = 2;

            face[0 + 1*(*nvf)] = 5;
            face[1 + 1*(*nvf)] = 6;
            face[2 + 1*(*nvf)] = 7;
            face[3 + 1*(*nvf)] = 8;

            face[0 + 2*(*nvf)] = 1;
            face[1 + 2*(*nvf)] = 2;
            face[2 + 2*(*nvf)] = 6;
            face[3 + 2*(*nvf)] = 5;

            face[0 + 3*(*nvf)] = 3;
            face[1 + 3*(*nvf)] = 4;
            face[2 + 3*(*nvf)] = 8;
            face[3 + 3*(*nvf)] = 7;

            face[0 + 4*(*nvf)] = 2;
            face[1 + 4*(*nvf)] = 3;
            face[2 + 4*(*nvf)] = 7;
            face[3 + 4*(*nvf)] = 6;

            face[0 + 5*(*nvf)] = 4;
            face[1 + 5*(*nvf)] = 1;
            face[2 + 5*(*nvf)] = 5;
            face[3 + 5*(*nvf)] = 8;
        }
    }
    else {
        fprintf(stderr, "Error: Only dim = 1, 2, or 3 supported.\n");
        exit(1);
    }
}

// Inputs:
// - int* t:       [nve * ne], element-to-node connectivity (column-major)
// - dstype* p:    [dim * np], coordinates of all points (column-major)
// - int* localfaces: [nvf * nfe], zero-based local face connectivity
// Outputs:
// - int* t2lf:    [nvf * nfe * ne], allocated by caller
// - dstype* pf:   [dim * nvf * nfe * ne], allocated by caller

void compute_pf(int* t2lf, dstype* pf,        
    const int* t, const dstype* p, const int* localfaces, 
    int nve, int nvf, int nfe, int dim, int ne, int np   
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
bool eval_expr(const char* expr_str, const dstype* x, int dim) {
    dstype x0 = x[0], x1 = (dim > 1 ? x[1] : 0), x2 = (dim > 2 ? x[2] : 0);
    te_variable vars[] = {
        {"x", &x0}, {"y", &x1}, {"z", &x2}
    };
    int err;
    te_expr* expr = te_compile(expr_str, vars, dim, &err);
    if (!expr) {
        fprintf(stderr, "TinyExpr error at %d in expression: %s\n", err, expr_str);
        return false;
    }
    dstype result = te_eval(expr);
    te_free(expr);
    return result != 0.0;
}

void assign_boundary_faces_with_expr(
    int* f,                     // [nfe x ne], zero-based, must be initialized to -1
    const int* f2t,             // [4 x nf], face-to-element connectivity
    const int* ind,             // [nb], indices of boundary faces   
    const dstype* pf,           // [dim x nvf x nfe x ne], column-major
    const char** bndexpr,       // boundary expressions, e.g., {"x*x + y*y <= 1"}
    int dim, int nvf, int nfe, int ne,
    int nb,                     // number of boundary faces        
    int nbndexpr
) {
    for (int i = 0; i < nfe * ne; ++i)
        f[i] = -1;

    int ind2 = (nvf > 1) ? (nvf - 1) : 1;

    for (int ii = 0; ii < nb; ++ii) {
        int b = ind[ii];
        int e = f2t[0 + 4 * b];
        int l = f2t[1 + 4 * b];

        for (int k = 0; k < nbndexpr; ++k) {
            bool ok = true;
            int check_pts[3] = {0, ind2, nvf - 1};

            for (int m = 0; m < 3; ++m) {
                int i = check_pts[m];
                const dstype* x = &pf[dim * (i + nvf * (l + nfe * e))];
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

void average_pf_over_face_points(
    dstype* pf,      // [dim * nvf * nfe * ne], column-major input/output
    int dim,
    int nvf,
    int nfe,
    int ne
) {
    dstype* temp = (dstype*)malloc(sizeof(dstype) * dim * nfe * ne);

    for (int e = 0; e < ne; ++e) {
        for (int f = 0; f < nfe; ++f) {
            for (int d = 0; d < dim; ++d) {
                dstype sum = 0.0;
                for (int v = 0; v < nvf; ++v) {
                    int idx = d + dim * (v + nvf * (f + nfe * e));
                    sum += pf[idx];
                }
                int out_idx = d + dim * (f + nfe * e);
                temp[out_idx] = sum / nvf;
            }
        }
    }

    // Copy result back to pf: now pf is shape [dim][nfe][ne]
    for (int i = 0; i < dim * nfe * ne; ++i)
        pf[i] = temp[i];

    free(temp);
}

void assign_periodic_bc_tinyexpr(    
    int* f, dstype* pf, const int* t, int* tprd, int* t2t,
    const int* t2lf,  const dstype* p,
    const int* prd_f1, const int* prd_f2,
    const char** expr1, const char** expr2,
    int nprd, int dim, int nfe, int ne, int nve, int nvf    
) {
    memcpy(tprd, t, sizeof(int) * nve * ne);

    for (int i = 0; i < nprd; ++i) {
        int* i1 = (int*)malloc(nfe * ne * sizeof(int));
        int* i2 = (int*)malloc(nfe * ne * sizeof(int));
        int n1 = 0, n2 = 0;

        for (int e = 0; e < ne; ++e) {
            for (int l = 0; l < nfe; ++l) {
                int idx = l + nfe * e;
                if (f[idx] == prd_f1[i]) {
                    f[idx] = -f[idx];
                    i1[n1++] = idx;
                } else if (f[idx] == prd_f2[i]) {
                    f[idx] = -f[idx];
                    i2[n2++] = idx;
                }
            }
        }

        dstype* p1 = (dstype*)malloc(dim * n1 * sizeof(dstype));
        dstype* p2 = (dstype*)malloc(dim * n2 * sizeof(dstype));

        for (int j = 0; j < n1; ++j)
            for (int d = 0; d < dim; ++d)
                p1[j * dim + d] = pf[d + dim * i1[j]];
        for (int j = 0; j < n2; ++j)
            for (int d = 0; d < dim; ++d)
                p2[j * dim + d] = pf[d + dim * i2[j]];

        dstype x, y, z;
        te_variable vars[] = {
            {"x", &x}, {"y", &y}, {"z", &z},
        };
        int len = (dim >= 3) ? 3 : dim;

        te_expr* ex1 = te_compile(expr1[i], vars, len, NULL);
        te_expr* ex2 = te_compile(expr2[i], vars, len, NULL);

        for (int j = 0; j < n1; ++j) {
            x = (dim > 0) ? p1[j * dim + 0] : 0;
            y = (dim > 1) ? p1[j * dim + 1] : 0;
            z = (dim > 2) ? p1[j * dim + 2] : 0;
            dstype val = te_eval(ex1);
            for (int d = 0; d < dim; ++d)
                p1[j * dim + d] = val; // Overwrite with evaluated (e.g., projected) result if needed
        }
        for (int j = 0; j < n2; ++j) {
            x = (dim > 0) ? p2[j * dim + 0] : 0;
            y = (dim > 1) ? p2[j * dim + 1] : 0;
            z = (dim > 2) ? p2[j * dim + 2] : 0;
            dstype val = te_eval(ex2);
            for (int d = 0; d < dim; ++d)
                p2[j * dim + d] = val;
        }

        te_free(ex1);
        te_free(ex2);

        int* e1 = (int*)malloc(n1 * sizeof(int));
        int* l1 = (int*)malloc(n1 * sizeof(int));
        int* e2 = (int*)malloc(n2 * sizeof(int));
        int* l2 = (int*)malloc(n2 * sizeof(int));

        for (int j = 0; j < n1; ++j) {
            e1[j] = i1[j] / nfe;
            l1[j] = i1[j] % nfe;
        }
        for (int j = 0; j < n2; ++j) {
            e2[j] = i2[j] / nfe;
            l2[j] = i2[j] % nfe;
        }

        int* in = (int*)malloc(n1 * sizeof(int));
        xiny<dstype>(in, p2, p1, n2, n1, dim);

        for (int j = 0; j < n1; ++j) {
            t2t[l2[in[j]] + nfe * e2[in[j]]] = e1[j];
            t2t[l1[j] + nfe * e1[j]] = e2[in[j]];
        }

        free(i1); free(i2); free(p1); free(p2);
        free(e1); free(l1); free(e2); free(l2); free(in);
    }
}

// 
// void assign_boundary_faces(
//     int* f,                  // Output: [nfe x ne], face markers (must be pre-allocated)
//     const int* f2t,          // [4 x nf], face-to-element connectivity
//     const int* ind,          // [nb], indices of boundary faces   
//     const dstype* pf,        // [dim x nvf x nfe x ne], column-major layout
//     int dim, int nvf, int nfe, int ne,
//     int nb,                  // number of boundary faces    
//     bool (**bndexpr)(const dstype*), // array of function pointers
//     int nbndexpr             // number of boundary expressions
// ) {
//     // Initialize f to -1
//     for (int i = 0; i < nfe * ne; ++i)
//         f[i] = -1;
// 
//     // Corrected interpretation of ind2 = min(1, nvf-1)
//     int ind2 = (nvf > 1) ? 1 : nvf - 1;
// 
//     for (int ii = 0; ii < nb; ++ii) {
//         int b = ind[ii];
//         int e = f2t[0 + 4 * b]; // element index
//         int l = f2t[1 + 4 * b]; // local face index
// 
//         for (int k = 0; k < nbndexpr; ++k) {
//             bool ok = true;
//             int check_pts[3] = {0, ind2, nvf - 1};
// 
//             for (int m = 0; m < 3; ++m) {
//                 int i = check_pts[m];
//                 const dstype* x = &pf[dim * (i + nvf * (l + nfe * e))];
//                 if (!bndexpr[k](x)) {
//                     ok = false;
//                     break;
//                 }
//             }
// 
//             if (ok) {
//                 f[l + nfe * e] = k;  // assign boundary marker
//                 break;
//             }
//         }
//     }
// }

// // Main function to assign periodic boundary conditions
// void assign_periodic_bc(
//     int nprd,
//     int dim,
//     int nfe,
//     int ne,
//     int nve,
//     int nvf,    
//     int* f,                 // [nfe * ne], face markers (modified)
//     dstype* pf,             // [dim * nfe * ne], face centers (modified)
//     const int* t,          // [nve * ne], element-to-node connectivity
//     int* tprd,             // [nve * ne], permuted connectivity (output)
//     int* t2t,              // [nfe * ne], face-to-element connectivity (modified)
//     const int* t2lf,       // [nvf * nface], face-to-node connectivity    
//     const dstype* p,       // [dim * np], node coordinates
//     const int* prd_f1,     // [nprd], face IDs for first periodic side
//     const int* prd_f2,     // [nprd], face IDs for second periodic side
//     void (*map1[])(const dstype*, int, dstype*), // prdexpr{i,2}
//     void (*map2[])(const dstype*, int, dstype*)  // prdexpr{i,4}
// ) {
//     memcpy(tprd, t, sizeof(int) * nve * ne);
// 
//     for (int i = 0; i < nprd; ++i) {
//         int* i1 = (int*)malloc(nfe * ne * sizeof(int));
//         int* i2 = (int*)malloc(nfe * ne * sizeof(int));
//         int n1 = 0, n2 = 0;
// 
//         for (int e = 0; e < ne; ++e) {
//             for (int l = 0; l < nfe; ++l) {
//                 int idx = l + nfe * e;
//                 if (f[idx] == prd_f1[i]) {
//                     f[idx] = -f[idx];
//                     i1[n1++] = idx;
//                 } else if (f[idx] == prd_f2[i]) {
//                     f[idx] = -f[idx];
//                     i2[n2++] = idx;
//                 }
//             }
//         }
// 
//         dstype* p1 = (dstype*)malloc(dim * n1 * sizeof(dstype));
//         dstype* p2 = (dstype*)malloc(dim * n2 * sizeof(dstype));
// 
//         for (int j = 0; j < n1; ++j)
//             for (int d = 0; d < dim; ++d)
//                 p1[j * dim + d] = pf[d + dim * i1[j]];
//         for (int j = 0; j < n2; ++j)
//             for (int d = 0; d < dim; ++d)
//                 p2[j * dim + d] = pf[d + dim * i2[j]];
// 
//         map1[i](p1, n1, p1);
//         map2[i](p2, n2, p2);
// 
//         int* e1 = (int*)malloc(n1 * sizeof(int));
//         int* l1 = (int*)malloc(n1 * sizeof(int));
//         int* e2 = (int*)malloc(n2 * sizeof(int));
//         int* l2 = (int*)malloc(n2 * sizeof(int));
// 
//         for (int j = 0; j < n1; ++j) {
//             e1[j] = i1[j] / nfe;
//             l1[j] = i1[j] % nfe;
//         }
//         for (int j = 0; j < n2; ++j) {
//             e2[j] = i2[j] / nfe;
//             l2[j] = i2[j] % nfe;
//         }
// 
//         int* in = (int*)malloc(n1 * sizeof(int));
//         xiny(p2, n2, p1, n1, dim, in);
// 
//         for (int j = 0; j < n1; ++j) {
//             t2t[l2[in[j]] + nfe * e2[in[j]]] = e1[j];
//             t2t[l1[j] + nfe * e1[j]] = e2[in[j]];
//         }
// 
//         int* v1 = (int*)malloc(nvf * n1 * sizeof(int));
//         int* v2 = (int*)malloc(nvf * n1 * sizeof(int));
//         for (int j = 0; j < n1; ++j)
//             for (int k = 0; k < nvf; ++k)
//                 v1[k + nvf * j] = t2lf[k + nvf * i1[j]];
//         for (int j = 0; j < n1; ++j)
//             for (int k = 0; k < nvf; ++k)
//                 v2[k + nvf * j] = t2lf[k + nvf * i2[in[j]]];
// 
//         int* uv1 = (int*)malloc(nve * sizeof(int));
//         int* uv2 = (int*)malloc(nve * sizeof(int));
//         memcpy(uv1, v1, sizeof(int) * nvf * n1);
//         memcpy(uv2, v2, sizeof(int) * nvf * n1);
//         int nv1 = unique_ints(uv1, nvf * n1);
//         int nv2 = unique_ints(uv2, nvf * n1);
// 
//         dstype* pp1 = (dstype*)malloc(nv1 * dim * sizeof(dstype));
//         dstype* pp2 = (dstype*)malloc(nv2 * dim * sizeof(dstype));
//         for (int j = 0; j < nv1; ++j)
//             for (int d = 0; d < dim; ++d)
//                 pp1[j * dim + d] = p[d + dim * uv1[j]];
//         for (int j = 0; j < nv2; ++j)
//             for (int d = 0; d < dim; ++d)
//                 pp2[j * dim + d] = p[d + dim * uv2[j]];
// 
//         map1[i](pp1, nv1, pp1);
//         map2[i](pp2, nv2, pp2);
// 
//         int* match = (int*)malloc(nv1 * sizeof(int));
//         xiny(pp2, nv2, pp1, nv1, dim, match);
// 
//         for (int j = 0; j < nv1; ++j)
//             for (int a = 0; a < ne * nve; ++a)
//                 if (tprd[a] == uv2[match[j]])
//                     tprd[a] = uv1[j];
// 
//         free(i1); free(i2); free(p1); free(p2);
//         free(e1); free(l1); free(e2); free(l2); free(in);
//         free(v1); free(v2); free(uv1); free(uv2);
//         free(pp1); free(pp2); free(match);
//     }
// }

// All arrays are assumed to be column-major, 0-based indexing
// p    : dstype[nd * npv]            -- flattened 2D array, p[dim + nd * j]
// t    : int[nve * ne]               -- connectivity: t[node + nve * e]
// philocal: dstype[npl * nve]        -- shape function values: phi[i + npl * node]
// dgnodes : dstype[npl * nd * ne]    -- output array
void compute_dgnodes(int npl, int nd, int ne, int nve, const dstype* philocal,
                     const dstype* p, const int* t, dstype* dgnodes)
{
    // Zero out dgnodes
    for (int i = 0; i < npl * nd * ne; ++i)
        dgnodes[i] = 0.0;

    for (int dim = 0; dim < nd; ++dim) {
        for (int node = 0; node < nve; ++node) {
            for (int e = 0; e < ne; ++e) {
                int global_node = t[node + nve * e];
                dstype pdim = p[dim + nd * global_node];  // p(dim, t(node,e))

                for (int i = 0; i < npl; ++i) {
                    dstype phi = philocal[i + npl * node];
                    int idx = i + npl * (dim + nd * e);   // dgnodes(i,dim,e)
                    dgnodes[idx] += phi * pdim;
                }
            }
        }
    }
}

void project_dgnodes_onto_curved_boundaries(
    int nd, int porder, int npe, int npf, int nfe, int ne,
    int* f, int* perm, int* curvedboundary,
    dstype* dgnodes,                         // shape: [npe][nd][ne] → column-major
    const char** fd_exprs                   // array of string expressions: fd_exprs[k] = "x*x + y*y - 1"
) {
    if (porder <= 1) return;

    int has_curved = 0;
    for (int k = 0; k < nfe * ne; ++k)
        if (f[k] > -1 && curvedboundary[f[k]] != 0)
            has_curved = 1;

    if (!has_curved) return;

    printf("Project dgnodes onto the curved boundaries...\n");

    dstype eps = 1e-14;
    dstype deps;

    for (int e = 0; e < ne; ++e) {
        for (int j = 0; j < nfe; ++j) {
            int fid = f[j + nfe * e];
            if (fid < 0) continue;

            int k = fid;
            if (curvedboundary[k] != 1) continue;

            const char* expr_str = fd_exprs[k];
            dstype x, y, z;
            te_variable vars[] = { {"x", &x}, {"y", &y}, {"z", &z} };
            int err;
            te_expr* expr = te_compile(expr_str, vars, nd, &err);
            if (!expr) {
                fprintf(stderr, "TinyExpr compile error at %d for expression: %s\n", err, expr_str);
                continue;
            }

            // collect coordinates: p[npts][nd]
            int npts = npf;
            dstype* p = (dstype*) malloc(sizeof(dstype) * npts * nd);
            for (int i = 0; i < npts; ++i) {
                int idx = perm[i + npf * j];
                for (int d = 0; d < nd; ++d)
                    p[i + npts * d] = dgnodes[idx + npe * (d + nd * e)];
            }

            // estimate deps
            dstype pmax = -1e20, pmin = 1e20;
            for (int i = 0; i < npts * nd; ++i) {
                if (p[i] > pmax) pmax = p[i];
                if (p[i] < pmin) pmin = p[i];
            }
            deps = sqrt(eps) * (pmax - pmin);

            // evaluate distance function d[i] = φ(p_i)
            dstype* d = (dstype*) malloc(sizeof(dstype) * npts);
            for (int i = 0; i < npts; ++i) {
                x = (nd > 0) ? p[i + npts * 0] : 0;
                y = (nd > 1) ? p[i + npts * 1] : 0;
                z = (nd > 2) ? p[i + npts * 2] : 0;
                d[i] = te_eval(expr);
            }

            // evaluate finite-difference gradients
            dstype* dgrad = (dstype*) malloc(sizeof(dstype) * npts * nd);
            for (int ddir = 0; ddir < nd; ++ddir) {
                for (int i = 0; i < npts; ++i)
                    p[i + npts * ddir] += deps;

                for (int i = 0; i < npts; ++i) {
                    x = (nd > 0) ? p[i + npts * 0] : 0;
                    y = (nd > 1) ? p[i + npts * 1] : 0;
                    z = (nd > 2) ? p[i + npts * 2] : 0;
                    dstype dshift = te_eval(expr);
                    dgrad[i + npts * ddir] = (dshift - d[i]) / deps;
                }

                for (int i = 0; i < npts; ++i)
                    p[i + npts * ddir] -= deps;
            }

            // project points onto the boundary
            for (int i = 0; i < npts; ++i) {
                dstype sumsq = 0;
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

            te_free(expr);
            free(d);
            free(dgrad);
            free(p);
        }
    }
}

int equal_row(const dstype* a, const dstype* b, int dim, dstype tol) {
    for (int i = 0; i < dim; ++i)
        if (fabs(a[i] - b[i]) > tol)
            return 0;
    return 1;
}

void get_node_coord(dstype* out, const dstype* dgnodes, int npe, int dim, int e, int a) {
    for (int d = 0; d < dim; ++d)
        out[d] = dgnodes[a + npe * (d + dim * e)];
}

/**
 * @param dgnodes   - [npe][dim][ne] in column-major layout
 * @param npe, dim, ne - mesh dimensions
 * @param cgnodes   - output array [max_nodes][dim], preallocated by user
 * @param cgelcon   - output array [npe][ne], preallocated by user
 * @return number of unique CG nodes (i.e., used portion of cgnodes)
 */
int mkelconcg(
    const dstype* dgnodes,  // [npe][dim][ne]
    int npe, int dim, int ne,
    dstype* cgnodes,        // [max_nodes][dim], user-allocated
    int* cgelcon            // [npe][ne], user-allocated
) {
    const dstype tol = 1e-12;
    int ncg = 0;  // number of unique CG nodes
    dstype node[3];  // supports up to 3D

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
                for (int d = 0; d < dim; ++d)
                    cgnodes[ncg * dim + d] = node[d];
                cgelcon[idx] = ncg;
                ncg++;
            }
        }
    }

    return ncg;  // return number of unique CG nodes
}

// Converts element-to-node connectivity (elcon) to node-to-element connectivity (CRS format)
// elcon: [nrow x ne] column-major, zero-based node indices
void mkent2elem(
    const int* elcon, // [nrow * ne], column-major
    int nrow,
    int ne,
    int** rowent2elem_out, // output [ndof+1], allocated inside
    int** colent2elem_out  // output [nnz], allocated inside
) {
    int total = nrow * ne;
    int entmax = 0;
    for (int i = 0; i < total; ++i)
        if (elcon[i] > entmax)
            entmax = elcon[i];

    bool* mark = (bool*)calloc(entmax + 1, sizeof(bool));
    for (int i = 0; i < total; ++i)
        if (elcon[i] >= 0)
            mark[elcon[i]] = true;

    int* ent2ind = (int*)malloc((entmax + 1) * sizeof(int));
    int ndof = 0;
    for (int i = 0; i <= entmax; ++i) {
        if (mark[i]) {
            ent2ind[i] = ndof++;
        } else {
            ent2ind[i] = -1;
        }
    }

    int* rowent2elem = (int*)calloc(ndof + 1, sizeof(int));
    int* counter = (int*)calloc(ndof, sizeof(int));
    bool* seen = (bool*)calloc(entmax + 1, sizeof(bool));

    for (int e = 0; e < ne; ++e) {
        for (int i = 0; i < nrow; ++i) {
            int ent = elcon[i + nrow * e];
            if (ent >= 0 && !seen[ent]) {
                seen[ent] = true;
                rowent2elem[ent2ind[ent] + 1]++;
            }
        }
        for (int i = 0; i < nrow; ++i) {
            int ent = elcon[i + nrow * e];
            if (ent >= 0)
                seen[ent] = false;
        }
    }

    for (int i = 1; i <= ndof; ++i)
        rowent2elem[i] += rowent2elem[i - 1];

    int* colent2elem = (int*)malloc(rowent2elem[ndof] * sizeof(int));

    for (int e = 0; e < ne; ++e) {
        for (int i = 0; i < nrow; ++i) {
            int ent = elcon[i + nrow * e];
            if (ent >= 0 && !seen[ent]) {
                seen[ent] = true;
                int id = ent2ind[ent];
                colent2elem[rowent2elem[id] + counter[id]++] = e;
            }
        }
        for (int i = 0; i < nrow; ++i) {
            int ent = elcon[i + nrow * e];
            if (ent >= 0)
                seen[ent] = false;
        }
    }

    free(mark);
    free(counter);
    free(seen);
    free(ent2ind);

    *rowent2elem_out = rowent2elem;
    *colent2elem_out = colent2elem;
}

// Helper function to match a point xcg to a row in xdg: returns 0-based index
int xiny(const dstype* xcg, const dstype* xdg, int npe, int dim) {
    for (int i = 0; i < npe; ++i) {
        int match = 1;
        for (int d = 0; d < dim; ++d) {
            if (fabs(xcg[d] - xdg[d + dim * i]) > 1e-12) {
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
    const int* rowent2elem, // [nent+1]
    const int* colent2elem, // [nnz]
    int nent,
    const dstype* cgnodes,  // [nent * dim], row-major
    const dstype* dgnodes,  // [npe * dim * ne], column-major
    int npe,
    int dim,
    int* cgent2dgent        // [nnz], output mapping (same shape as colent2elem)
) {
    for (int i = 0; i < nent; ++i) {
        int start = rowent2elem[i];
        int end = rowent2elem[i + 1];
        const dstype* xcg = &cgnodes[i * dim];
        for (int j = start; j < end; ++j) {
            int elem = colent2elem[j];
            const dstype* xdg = &dgnodes[npe * dim * elem];
            int in = xiny(xcg, xdg, npe, dim);
            cgent2dgent[j] = elem * npe + in; // global DG index
        }
    }
}
