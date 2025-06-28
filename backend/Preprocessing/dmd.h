#include <vector>
#include <algorithm>
#include <numeric>

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

// Helper to match rows of A and B with tolerance; returns indices (1-based)
void xiny(const double* A, const double* B, int npf, int dim, int* result) 
{
    int tol = 1e-8;
    for (int i = 0; i < npf; ++i) {
        int found = 0;
        for (int j = 0; j < npf; ++j) {
            int match = 1;
            for (int d = 0; d < dim; ++d) {
                if (fabs(B[i + d*npf] - A[j + d*npf]) > tol) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                result[i] = j + 1;  // 1-based
                found = 1;
                break;
            }
        }
        if (!found) result[i] = 0;  // not found
    }
}

// Main permindex function
int* permindex(const int* plocfc, int npf, int dim, int elemtype, int* ncols_out) {
    int* ind = NULL;

    if (dim == 1) {
        *ncols_out = 1;
        ind = (int*) malloc(sizeof(int));
        ind[0] = 1;
    } 
    else if (dim == 2) {
        *ncols_out = 1;
        ind = (int*) malloc(sizeof(int) * npf);
        for (int i = 0; i < npf; ++i)
            ind[i] = npf - i;
    } 
    else if (dim == 3 && elemtype == 0) {
        *ncols_out = 3;
        ind = (int*) malloc(sizeof(int) * npf * 3);
        int* plocfc2 = (int*) malloc(sizeof(int) * npf * 2);

        // Perm 1: swap columns
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 0*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 0*npf);

        // Perm 2: 1 - xi - eta in col 1
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 0*npf] - plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 1*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 1*npf);

        // Perm 3: 1 - xi - eta in col 2
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 0*npf] - plocfc[i + 1*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 2*npf);

        free(plocfc2);
    } 
    else if (dim == 3 && elemtype == 1) {
        *ncols_out = 4;
        ind = (int*) malloc(sizeof(int) * npf * 4);
        int* plocfc2 = (int*) malloc(sizeof(int) * npf * 2);

        // Perm 1: swap columns
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = plocfc[i + 0*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 0*npf);

        // Perm 2: eta = 1 - eta
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 1*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 1*npf);

        // Perm 3: xi = 1 - eta, eta = 1 - xi
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 1*npf];
            plocfc2[i + 1*npf] = 1.0 - plocfc[i + 0*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 2*npf);

        // Perm 4: xi = 1 - xi
        for (int i = 0; i < npf; ++i) {
            plocfc2[i + 0*npf] = 1.0 - plocfc[i + 0*npf];
            plocfc2[i + 1*npf] = plocfc[i + 1*npf];
        }
        xiny(plocfc, plocfc2, npf, 2, ind + 3*npf);

        free(plocfc2);
    }

    return ind;
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

            if (e2 > 0) {
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
        for (int i = 0; i < nf; ++i) {
            int e1 = f2t[0 + 4 * i];
            int l1 = f2t[1 + 4 * i] - 1;
            int e2 = f2t[2 + 4 * i];
            int l2 = f2t[3 + 4 * i] - 1;

            for (int j = 0; j < npf; ++j)
                facecon[0 + 2 * j + 2 * npf * i] = e1 * npe + perm[j + npf * l1];

            for (int j = 0; j < npf; ++j)
                elemcon[j + npf * (l1 + nfe * (e1 - 1))] = i * npf + facenode1[j];

            if (e2 > 0) {
                int f1[4], f2[4];
                for (int k = 0; k < (elemtype == 0 ? 3 : 4); ++k) {
                    f1[k] = t[face[k + 4 * l1] + ne * e1];
                    f2[k] = t[face[k + 4 * l2] + ne * e2];
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
                    fprintf(stderr, "Mesh connectivity is wrong\n");
                    exit(1);
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

    free(facenode1);
    free(facenode2);
}

void apply_bcm(int* bf, const int* fi, const int* bcm, size_t n, size_t nbcm) 
{
    // initialize bf to 0
    for (size_t i = 0; i < n; ++i)
        bf[i] = 0;

    // apply boundary condition map
    for (size_t j = 0; j < nbcm; ++j) {
        for (size_t i = 0; i < n; ++i) {
            if (fi[i] == (int)(j)) {  // MATLAB is 1-based
                bf[i] = bcm[j];
            }
        }
    }
}

// Copy selected columns of a column-major matrix
// a: input matrix (m × n), column-major
// ind: list of column indices to keep (length k)
// a_new: output matrix (m × k), pre-allocated
void select_columns(int* a_new, const int* a, const int* ind, int m, int n, int k) 
{
    for (int j = 0; j < k; ++j) {
        int col = ind[j];
        for (int i = 0; i < m; ++i) {
            a_new[i + j * m] = a[i + col * m];
        }
    }
}

void extract_fb(int* fb, const int* inb, const int* f2t, const int* fi, int nb, int nfe) 
{
    for (int i = 0; i < nb; ++i) {
        int b  = inb[i];                  // MATLAB: b = inb(i)
        int e1 = f2t[0 + 4 * b];          // MATLAB: f2t(1,b)
        int l1 = f2t[1 + 4 * b];          // MATLAB: f2t(2,b)
        fb[i] = fi[l1 + nfe * e1];        // MATLAB: fi(l1,e1)
    }
}

void extract_subset(int* b, const int* a, const int* ind, size_t k) 
{
    for (size_t i = 0; i < k; ++i) {
        b[i] = a[ind[i]];
    }
}

void sort_with_indices(int* b, int* ind, const int* a, int n) 
{
    // Initialize b and ind
    for (int i = 0; i < n; ++i) {
        b[i] = a[i];
        ind[i] = i;
    }

    // Selection sort with index tracking
    for (int i = 0; i < n - 1; ++i) {
        int min_idx = i;
        for (int j = i + 1; j < n; ++j) {
            if (b[j] < b[min_idx]) {
                min_idx = j;
            }
        }
        // Swap b[i] and b[min_idx]
        if (min_idx != i) {
            double tmp_val = b[i];
            b[i] = b[min_idx];
            b[min_idx] = tmp_val;

            int tmp_idx = ind[i];
            ind[i] = ind[min_idx];
            ind[min_idx] = tmp_idx;
        }
    }
}

#include <stddef.h>

/*-----------------------------------------------------------
 * unique_count
 *  a   : input sorted array of length n   (double)
 *  n   : length of a
 *  b   : output array (unique values)     (double)
 *  c   : output array (occurrence counts) (int)
 *  return value = number of unique entries written to b & c
 *
 *  Caller must pre-allocate b and c with size >= n.
 *----------------------------------------------------------*/
size_t unique_count(int *b, int *c, const int *a, size_t n)
{
    if (n == 0) return 0;

    size_t uniq = 0;          /* index in b/c */
    double current = a[0];
    int     count   = 1;

    for (size_t i = 1; i < n; ++i) {
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

void fix_f2t_consistency(int* f2t, const int* elempart, int nf) {
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

void localbasis(int *phielem, int *phiface, const int *plocvl, const int *plocfc,
                int dim, int elemtype, int nne, int nnf)
{
    int i;

    if (dim == 1) {  // 1D line element
        for (i = 0; i < nnf; ++i) {
            phiface[i] = 1.0;  // scalar constant
        }

        for (i = 0; i < nne; ++i) {
            int xi = plocvl[i];  // since plocvl(:,1)
            phielem[i + 0 * nne] = 1.0 - xi;  // phielem(:,1)
            phielem[i + 1 * nne] = xi;        // phielem(:,2)
        }
    }
    else if (dim == 2 && elemtype == 0) {  // triangle
        for (i = 0; i < nnf; ++i) {
            int xi = plocfc[i];  // plocfc(:,1)
            phiface[i + 0 * nnf] = 1.0 - xi;  // phiface(:,1)
            phiface[i + 1 * nnf] = xi;        // phiface(:,2)
        }

        for (i = 0; i < nne; ++i) {
            int xi  = plocvl[i + 0 * nne];  // plocvl(:,1)
            int eta = plocvl[i + 1 * nne];  // plocvl(:,2)
            phielem[i + 0 * nne] = 1.0 - xi - eta;
            phielem[i + 1 * nne] = xi;
            phielem[i + 2 * nne] = eta;
        }
    }
    else if (dim == 2 && elemtype == 1) {  // quadrilateral
        for (i = 0; i < nnf; ++i) {
            int xi = plocfc[i + 0 * nnf];
            phiface[i + 0 * nnf] = 1.0 - xi;
            phiface[i + 1 * nnf] = xi;
        }

        for (i = 0; i < nne; ++i) {
            int xi  = plocvl[i + 0 * nne];
            int eta = plocvl[i + 1 * nne];
            phielem[i + 0 * nne] = (1.0 - xi) * (1.0 - eta);
            phielem[i + 1 * nne] =  xi        * (1.0 - eta);
            phielem[i + 2 * nne] =  xi        * eta;
            phielem[i + 3 * nne] = (1.0 - xi) * eta;
        }
    }
    else if (dim == 3 && elemtype == 0) {  // tetrahedron
        for (i = 0; i < nnf; ++i) {
            int xi  = plocfc[i + 0 * nnf];
            int eta = plocfc[i + 1 * nnf];
            phiface[i + 0 * nnf] = 1.0 - xi - eta;
            phiface[i + 1 * nnf] = xi;
            phiface[i + 2 * nnf] = eta;
        }

        for (i = 0; i < nne; ++i) {
            int xi   = plocvl[i + 0 * nne];
            int eta  = plocvl[i + 1 * nne];
            int zeta = plocvl[i + 2 * nne];
            phielem[i + 0 * nne] = 1.0 - xi - eta - zeta;
            phielem[i + 1 * nne] = xi;
            phielem[i + 2 * nne] = eta;
            phielem[i + 3 * nne] = zeta;
        }
    }
    else if (dim == 3 && elemtype == 1) {  // hexahedron
        for (i = 0; i < nnf; ++i) {
            int xi  = plocfc[i + 0 * nnf];
            int eta = plocfc[i + 1 * nnf];
            phiface[i + 0 * nnf] = (1.0 - xi) * (1.0 - eta);
            phiface[i + 1 * nnf] =  xi        * (1.0 - eta);
            phiface[i + 2 * nnf] =  xi        * eta;
            phiface[i + 3 * nnf] = (1.0 - xi) * eta;
        }

        for (i = 0; i < nne; ++i) {
            int xi   = plocvl[i + 0 * nne];
            int eta  = plocvl[i + 1 * nne];
            int zeta = plocvl[i + 2 * nne];
            phielem[i + 0 * nne] = (1.0 - xi) * (1.0 - eta) * (1.0 - zeta);
            phielem[i + 1 * nne] =  xi        * (1.0 - eta) * (1.0 - zeta);
            phielem[i + 2 * nne] =  xi        * eta         * (1.0 - zeta);
            phielem[i + 3 * nne] = (1.0 - xi) * eta         * (1.0 - zeta);
            phielem[i + 4 * nne] = (1.0 - xi) * (1.0 - eta) * zeta;
            phielem[i + 5 * nne] =  xi        * (1.0 - eta) * zeta;
            phielem[i + 6 * nne] =  xi        * eta         * zeta;
            phielem[i + 7 * nne] = (1.0 - xi) * eta         * zeta;
        }
    }
}


void xiny(std::vector<int>& indices, const std::vector<int>& x, const std::vector<int>& y) {
    indices.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        auto it = std::find(y.begin(), y.end(), x[i]);
        indices[i] = (it != y.end()) ? std::distance(y.begin(), it) : -1;
    }
}

std::vector<int> neighboringelements(const int* e2e, const int* elem, int nfe, int ne, int nelem) {
    std::vector<int> neighbors;
    neighbors.reserve(nelem * nfe);

    for (int k = 0; k < nelem; ++k) {
        int e = elem[k];
        if (e < 0 || e >= ne) continue;
        for (int i = 0; i < nfe; ++i) {
            int nb = e2e[i + e * nfe];
            if (nb >= 0) {
                neighbors.push_back(nb);
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

    return neighbors;
}

void create_elemsend(std::vector<dmdstruct>& dmd);

void build_dmd(std::vector<dmdstruct>& dmd, const int* e2e, const int* elem2cpu, const int* inte, int nfe, int ne, int coupledinterface) {
    int nproc = static_cast<int>(dmd.size());

    for (int i = 0; i < nproc; ++i) {
        std::vector<int> intelem;
        intelem.reserve(ne / nproc);
        for (int e = 0; e < ne; ++e) {
            if (elem2cpu[e] == i) intelem.push_back(e);
        }

        std::vector<int> elem = neighboringelements(e2e, intelem.data(), nfe, ne, intelem.size());
        std::vector<int> extelem;
        extelem.reserve(elem.size());
        std::set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(extelem));

        elem = neighboringelements(e2e, extelem.data(), nfe, ne, extelem.size());
        std::vector<int> bndelem;
        bndelem.reserve(elem.size());
        std::set_intersection(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(bndelem));

        std::vector<int> bndelem1, bndelem2;
        if (coupledinterface > 0) {
            bndelem1.reserve(bndelem.size());
            bndelem2.reserve(bndelem.size());
            std::set_intersection(bndelem.begin(), bndelem.end(), inte, inte + ne, std::back_inserter(bndelem1));
            std::set_difference(bndelem.begin(), bndelem.end(), bndelem1.begin(), bndelem1.end(), std::back_inserter(bndelem2));
            bndelem = bndelem1;
            bndelem.insert(bndelem.end(), bndelem2.begin(), bndelem2.end());
        }

        std::vector<int> part1, part2;
        part1.reserve(intelem.size());
        part2 = bndelem;
        std::set_difference(intelem.begin(), intelem.end(), bndelem.begin(), bndelem.end(), std::back_inserter(part1));
        std::vector<int>& part3 = extelem;

        dmd[i].elempart.clear();
        dmd[i].elempart.reserve(part1.size() + part2.size() + part3.size());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part1.begin(), part1.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part2.begin(), part2.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part3.begin(), part3.end());

        dmd[i].elempartpts = {static_cast<int>(part1.size()), static_cast<int>(part2.size()), static_cast<int>(part3.size())};

        if (coupledinterface > 0) {
            dmd[i].intepartpts = {static_cast<int>(part1.size()), static_cast<int>(bndelem1.size()), static_cast<int>(bndelem2.size()), static_cast<int>(part3.size())};
        }

        dmd[i].elem2cpu.resize(dmd[i].elempart.size());
        for (size_t j = 0; j < dmd[i].elempart.size(); ++j)
            dmd[i].elem2cpu[j] = elem2cpu[dmd[i].elempart[j]];

        std::vector<int> recvelem = extelem;
        std::vector<int> ind;
        xiny(ind, recvelem, dmd[i].elempart);

        dmd[i].elemrecv.clear();
        dmd[i].elemrecv.reserve(recvelem.size());
        int offset = part1.size() + part2.size();
        for (size_t j = 0; j < recvelem.size(); ++j) {
            if (ind[j] >= 0) {
                dmd[i].elemrecv.push_back({dmd[i].elem2cpu[ind[j]] + 1, static_cast<int>(offset + j), recvelem[j]});
            }
        }

        std::sort(dmd[i].elemrecv.begin(), dmd[i].elemrecv.end());

        dmd[i].nbsd.clear();
        dmd[i].nbsd.reserve(dmd[i].elemrecv.size());
        for (const auto& row : dmd[i].elemrecv)
            dmd[i].nbsd.push_back(row[0]);

        std::sort(dmd[i].nbsd.begin(), dmd[i].nbsd.end());
        dmd[i].nbsd.erase(std::unique(dmd[i].nbsd.begin(), dmd[i].nbsd.end()), dmd[i].nbsd.end());
    }

    create_elemsend(dmd);
}

void build_dmd(std::vector<dmdstruct>& dmd, const int* e2e, const int* elem2cpu, int nfe, int ne) {
    int nproc = static_cast<int>(dmd.size());

    for (int i = 0; i < nproc; ++i) {
        std::vector<int> intelem;
        intelem.reserve(ne / nproc);
        for (int e = 0; e < ne; ++e) {
            if (elem2cpu[e] == i) intelem.push_back(e);
        }

        std::vector<int> elem = neighboringelements(e2e, intelem.data(), nfe, ne, intelem.size());

        std::vector<int> extelem;
        extelem.reserve(elem.size());
        std::set_difference(elem.begin(), elem.end(), intelem.begin(), intelem.end(), std::back_inserter(extelem));

        std::vector<int> elem2 = neighboringelements(e2e, extelem.data(), nfe, ne, extelem.size());

        std::vector<int> bndelem;
        bndelem.reserve(elem2.size());
        std::set_intersection(elem2.begin(), elem2.end(), intelem.begin(), intelem.end(), std::back_inserter(bndelem));

        std::vector<int> intex;
        std::set_union(intelem.begin(), intelem.end(), extelem.begin(), extelem.end(), std::back_inserter(intex));
        std::vector<int> outelem;
        outelem.reserve(elem2.size());
        std::set_difference(elem2.begin(), elem2.end(), intex.begin(), intex.end(), std::back_inserter(outelem));

        std::vector<int> part1, part2, part3, part4;
        part1.reserve(intelem.size());
        part2 = bndelem;
        part3 = extelem;
        part4 = outelem;

        std::set_difference(intelem.begin(), intelem.end(), bndelem.begin(), bndelem.end(), std::back_inserter(part1));

        dmd[i].elempart.clear();
        dmd[i].elempart.reserve(part1.size() + part2.size() + part3.size() + part4.size());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part1.begin(), part1.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part2.begin(), part2.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part3.begin(), part3.end());
        dmd[i].elempart.insert(dmd[i].elempart.end(), part4.begin(), part4.end());

        dmd[i].elempartpts = {static_cast<int>(part1.size()), static_cast<int>(part2.size()), static_cast<int>(part3.size()), static_cast<int>(part4.size())};

        dmd[i].elem2cpu.resize(dmd[i].elempart.size());
        for (size_t j = 0; j < dmd[i].elempart.size(); ++j)
            dmd[i].elem2cpu[j] = elem2cpu[dmd[i].elempart[j]];

        std::vector<int> recvelem = extelem;
        recvelem.insert(recvelem.end(), outelem.begin(), outelem.end());

        std::vector<int> ind;
        xiny(ind, recvelem, dmd[i].elempart);

        dmd[i].elemrecv.clear();
        dmd[i].elemrecv.reserve(recvelem.size());
        int offset = part1.size() + part2.size();
        for (size_t j = 0; j < recvelem.size(); ++j) {
            if (ind[j] >= 0) {
                dmd[i].elemrecv.push_back({dmd[i].elem2cpu[ind[j]] + 1, static_cast<int>(offset + j), recvelem[j]});
            }
        }

        std::sort(dmd[i].elemrecv.begin(), dmd[i].elemrecv.end());

        dmd[i].nbsd.clear();
        dmd[i].nbsd.reserve(dmd[i].elemrecv.size());
        for (const auto& row : dmd[i].elemrecv)
            dmd[i].nbsd.push_back(row[0]);

        std::sort(dmd[i].nbsd.begin(), dmd[i].nbsd.end());
        dmd[i].nbsd.erase(std::unique(dmd[i].nbsd.begin(), dmd[i].nbsd.end()), dmd[i].nbsd.end());
    }

    create_elemsend(dmd);
}

void create_elemsend(std::vector<dmdstruct>& dmd) {
    int nproc = dmd.size();

    for (int k = 0; k < nproc; ++k) {
        dmd[k].elemsend.clear();
    }

    for (int i = 0; i < nproc; ++i) {
        for (size_t j = 0; j < dmd[i].nbsd.size(); ++j) {
            int k = dmd[i].nbsd[j];

            std::vector<std::vector<int>> tm;
            tm.reserve(dmd[i].elemrecv.size());
            for (const auto& row : dmd[i].elemrecv) {
                if (row[0] == k) tm.push_back(row);
            }

            for (auto& row : tm) row[0] = i;

            std::vector<int> x(tm.size());
            for (size_t m = 0; m < tm.size(); ++m)
                x[m] = tm[m][2];

            std::vector<int> local_indices;
            local_indices.reserve(x.size());
            xiny(local_indices, x, dmd[k].elempart);

            for (size_t m = 0; m < tm.size(); ++m)
                tm[m][2] = local_indices[m];

            dmd[k].elemsend.insert(dmd[k].elemsend.end(), tm.begin(), tm.end());
        }
    }

    for (int i = 0; i < nproc; ++i) {
        dmd[i].elemsendpts.resize(dmd[i].nbsd.size());
        dmd[i].elemrecvpts.resize(dmd[i].nbsd.size());

        for (size_t j = 0; j < dmd[i].nbsd.size(); ++j) {
            int n = dmd[i].nbsd[j];

            dmd[i].elemsendpts[j] = std::count_if(
                dmd[i].elemsend.begin(), dmd[i].elemsend.end(),
                [n](const std::vector<int>& row) { return row[0] == n; });

            dmd[i].elemrecvpts[j] = std::count_if(
                dmd[i].elemrecv.begin(), dmd[i].elemrecv.end(),
                [n](const std::vector<int>& row) { return row[0] == n; });
        }
    }
}
