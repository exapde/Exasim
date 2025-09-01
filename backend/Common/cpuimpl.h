/**
 * @file cpuimpl.h
 * @brief Common CPU implementations for mesh and array operations in Exasim backend.
 *
 * This header provides a collection of utility functions for manipulating arrays,
 * mesh connectivity, face and element indexing, sorting, and CRS (Compressed Row Storage)
 * matrix construction. It is designed for finite element and finite volume mesh operations,
 * including face-to-element, element-to-face, and neighbor relationships.
 *
 * Main functionalities:
 * - Array operations: setting values, finding min/max, inserting subarrays.
 * - Mesh connectivity: building face-to-element (f2e), element-to-element (e2e),
 *   element-to-face (e2f), and face-to-face (f2f) mappings.
 * - Face and element indexing: generating indices for faces and elements, including
 *   boundary and interface faces.
 * - Sorting and ordering: bubble sort, recurrence-based sort, and reordering of faces
 *   for efficient computation.
 * - CRS matrix construction: building row pointers and column indices for sparse matrices
 *   representing mesh connectivity.
 * - Partitioning: grid partitioning for 2D and 3D meshes into blocks.
 * - ILU0 indexing: generating indices for incomplete LU factorization on CRS matrices.
 *
 * Types:
 * - dstype: Data type for array operations (typically double or float).
 * - Int: Integer type for printing arrays (typically int).
 *
 * Usage:
 * - Include this header in CPU-based mesh processing modules.
 * - Functions are designed for 0-based indexing and column-major storage where applicable.
 * - Memory allocation is performed internally for temporary arrays; caller is responsible
 *   for freeing output arrays where noted.
 *
 * Note:
 * - Many functions assume specific mesh layouts and are tailored for Exasim's mesh data structures.
 * - Error checking is performed for grid partitioning and ILU0 diagonal presence.
 *
 * Author: Exasim backend team
 * License: MIT or project-specific
 */
#ifndef __CPUIMPL_H__
#define __CPUIMPL_H__

void cpuArraySetValue(dstype *udg, dstype a, int N)
{
  for (int i=0; i<N; i++)
    udg[i] = a;
}

void cpuArrayInsert(dstype* u, const dstype* un, const int I, const int J, const int K, 
        const int i1, const int i2, const int j1, const int j2, const int k1, const int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int Q = I*J;    
    for (int idx=0; idx<N; idx++) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+Q*k] = un[idx];        
    }
}


dstype cpuArrayMin(dstype *a, int n)
{
    dstype b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

dstype cpuArrayMax(dstype *a, int n)
{
    dstype b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

void cpuElemGeom(dstype *Xx, dstype *jac, dstype *Jg, int ne, int ng, int nd)
{        
    dstype *Jg11, *Jg12, *Jg13, *Jg21, *Jg22, *Jg23, *Jg31, *Jg32, *Jg33;
    dstype *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;
    int ngv = ng*ne;
     
    if (nd==1) {
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg[i];
            Xx[i] = 1.0;
        }
    }
    else if (nd==2) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg21 = &Jg[2*ngv];
        Jg22 = &Jg[3*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx12 = &Xx[2*ngv];
        Xx22 = &Xx[3*ngv];
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
            Xx11[i] = Jg22[i]; // dxi/dx
            Xx21[i] = -Jg21[i]; //dxi/dy 
            Xx12[i] = -Jg12[i]; //deta/dx
            Xx22[i] = Jg11[i]; //deta/dy
        }        
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg12 = &Jg[ngv];
        Jg13 = &Jg[2*ngv];
        Jg21 = &Jg[3*ngv];
        Jg22 = &Jg[4*ngv];
        Jg23 = &Jg[5*ngv];
        Jg31 = &Jg[6*ngv];
        Jg32 = &Jg[7*ngv];
        Jg33 = &Jg[8*ngv];
        Xx11 = &Xx[0];
        Xx21 = &Xx[ngv];
        Xx31 = &Xx[2*ngv];
        Xx12 = &Xx[3*ngv];
        Xx22 = &Xx[4*ngv];
        Xx32 = &Xx[5*ngv];
        Xx13 = &Xx[6*ngv];
        Xx23 = &Xx[7*ngv];
        Xx33 = &Xx[8*ngv];
        for (int i=0; i<ngv; i++) {
            jac[i] = Jg11[i]*Jg22[i]*Jg33[i] - Jg11[i]*Jg32[i]*Jg23[i] +
                     Jg21[i]*Jg32[i]*Jg13[i] - Jg21[i]*Jg12[i]*Jg33[i] +
                     Jg31[i]*Jg12[i]*Jg23[i] - Jg31[i]*Jg22[i]*Jg13[i];
            Xx11[i] = Jg22[i]*Jg33[i] - Jg23[i]*Jg32[i];
            Xx21[i] = Jg23[i]*Jg31[i] - Jg21[i]*Jg33[i];
            Xx31[i] = Jg21[i]*Jg32[i] - Jg22[i]*Jg31[i];
            Xx12[i] = Jg13[i]*Jg32[i] - Jg12[i]*Jg33[i];
            Xx22[i] = Jg11[i]*Jg33[i] - Jg13[i]*Jg31[i];
            Xx32[i] = Jg12[i]*Jg31[i] - Jg11[i]*Jg32[i];
            Xx13[i] = Jg12[i]*Jg23[i] - Jg13[i]*Jg22[i];
            Xx23[i] = Jg13[i]*Jg21[i] - Jg11[i]*Jg23[i];
            Xx33[i] = Jg11[i]*Jg22[i] - Jg12[i]*Jg21[i];
        }        
    }
}

void cpuGetElemNodes(dstype* unView, const dstype* uView, const int np, const int nc, const int nc1, const int nc2, const int e1, const int e2) 
{
    int nn = np * (e2 - e1);
    int ncu = nc2 - nc1;
    int N = nn * ncu;
    int K = np * nc;

    for (int idx=0; idx<N; idx++) {
        int i = idx % nn;  // [0, np*(e2-e1)]
        int j = idx / nn;  // [0, ncu]
        int k = i % np;    // [0, np]
        int e = i / np + e1;
        unView[idx] = uView[k + (j + nc1) * np + e * K];
    }
}

void cpuApplyGivensRotation(dstype *H, dstype *s, dstype *cs, dstype *sn,  int i)
{        
    dstype temp;    
    for (int k=0; k<i; k++) {
        temp       =  cs[k]*H[k] + sn[k]*H[k+1];
        H[k+1] = -sn[k]*H[k] + cs[k]*H[k+1];
        H[k]   = temp;
    }

    if (H[i+1] == 0.0) {
        cs[i] = 1.0;
        sn[i] = 0.0;
    } 
    else if (fabs(H[i+1]) > fabs(H[i])) {
        temp = H[i] / H[i+1];
        sn[i] = 1.0 / sqrt( 1.0 + temp*temp );
        cs[i] = temp * sn[i];
    } 
    else {
        temp = H[i+1] / H[i];
        cs[i] = 1.0 / sqrt( 1.0 + temp*temp );
        sn[i] = temp * cs[i];
    }
   
    temp   = cs[i]*s[i];                       
    s[i+1] = -sn[i]*s[i];
    s[i]   = temp;
    H[i] = cs[i]*H[i] + sn[i]*H[i+1];
    H[i+1] = 0.0;
}

void cpuBackSolve(dstype *y, dstype *H, dstype *s, int i, int n)
{
    for (int j=i; j>=0; j--)
        y[j] = s[j];    
    
    for (int j=i; j>=0; j--) {
        y[j] =  y[j]/H[j+n*j]; 
        for (int k=j-1; k>=0; k--)
            y[k] = y[k] - H[k+n*j]*y[j]; 
    }
}

int find(const int* a, int b, int n, int opts) 
{
    int count = 0;
    if (opts==0) {
      for (int i = 0; i < n; ++i) {
          if (a[i] == b) count++;                        
      }
    }
    else if (opts==1) {
      for (int i = 0; i < n; ++i) {
          if (a[i] <= b) count++;
      }
    }
    else if (opts==2) {
      for (int i = 0; i < n; ++i) {
          if (a[i] >= b) count++;
      }
    }    
    return count;
}

int find(int* indices, const int* a, int b, int n, int opts) 
{
    int count = 0;
    if (opts==0) {
      for (int i = 0; i < n; ++i) {
          if (a[i] == b) {
              indices[count++] = i;
          }
      }
    }
    else if (opts==1) {
      for (int i = 0; i < n; ++i) {
          if (a[i] <= b) {
              indices[count++] = i;
          }
      }
    }
    else if (opts==2) {
      for (int i = 0; i < n; ++i) {
          if (a[i] >= b) {
              indices[count++] = i;
          }
      }
    }
    
    return count;
}

bool equal_face_nodes(const int* a, const int* b, int nnf) {
    for (int i = 0; i < nnf; i++) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

// Sort small arrays (nnf elements) using bubble sort
void sort_face_nodes(int* face, int nnf) {
    for (int i = 0; i < nnf - 1; ++i) {
        for (int j = 0; j < nnf - i - 1; ++j) {
            if (face[j] > face[j + 1]) {
                int tmp = face[j];
                face[j] = face[j + 1];
                face[j + 1] = tmp;
            }
        }
    }
}

int mkf2e(int* f2e, const int* e2n, const int* local_faces, int ne, int nne, int nnf, int nfe) 
{
    int max_faces = ne * nfe;
    int* face_nodes = (int*)malloc(sizeof(int) * nnf * max_faces); // stores sorted face nodes
    int count = 0;

    for (int elem = 0; elem < ne; ++elem) {
        for (int lf = 0; lf < nfe; ++lf) {
            int temp_face[8];  
            for (int i = 0; i < nnf; ++i) {
                int local_node = local_faces[i + lf * nnf];
                temp_face[i] = e2n[elem * nne + local_node];
            }

            sort_face_nodes(temp_face, nnf);

            // Check if face already exists
            int match = -1;
            for (int f = 0; f < count; ++f) {
                if (equal_face_nodes(temp_face, &face_nodes[f * nnf], nnf)) {
                    match = f;
                    break;
                }
            }

            if (match >= 0) {
                // Second element sharing the face
                f2e[2 + match * 4] = elem;  // element on other side
                f2e[3 + match * 4] = lf;    // local face index in that element
            } else {
                // New face
                for (int i = 0; i < nnf; ++i)
                    face_nodes[count * nnf + i] = temp_face[i];

                f2e[0 + count * 4] = elem;  // first element
                f2e[1 + count * 4] = lf;    // local face index in first element
                f2e[2 + count * 4] = -1;    // no neighbor yet
                f2e[3 + count * 4] = -1;
                ++count;
            }
        }
    }
    
    free(face_nodes);
    return count;    
}

int mke2e(int* e2e, const int* e2n, const int* local_faces, int ne, int nne, int nnf, int nfe) 
{
    int max_faces = ne * nfe;
    int* face_nodes = (int*)malloc(sizeof(int) * nnf * max_faces);
    int* face_owner = (int*)malloc(sizeof(int) * max_faces);
    int* face_lf = (int*)malloc(sizeof(int) * max_faces);
    int face_count = 0;

    for (int elem = 0; elem < ne; ++elem) {
        for (int lf = 0; lf < nfe; ++lf) {
            int temp_face[8];
            for (int i = 0; i < nnf; ++i) {
                int local_node = local_faces[i + lf * nnf];
                temp_face[i] = e2n[elem * nne + local_node];
            }

            sort_face_nodes(temp_face, nnf);

            int match = -1;
            for (int f = 0; f < face_count; ++f) {
                if (equal_face_nodes(temp_face, &face_nodes[f * nnf], nnf)) {
                    match = f;
                    break;
                }
            }

            if (match >= 0) {
                // Face shared with another element
                int e1 = face_owner[match];
                int l1 = face_lf[match];
                int e2 = elem;
                int l2 = lf;

                e2e[e1 * nfe + l1] = e2;
                e2e[e2 * nfe + l2] = e1;
            } else {
                // New face
                for (int i = 0; i < nnf; ++i)
                    face_nodes[face_count * nnf + i] = temp_face[i];
                face_owner[face_count] = elem;
                face_lf[face_count] = lf;
                e2e[elem * nfe + lf] = -1;  // boundary face
                face_count++;
            }
        }
    }

    free(face_nodes);
    free(face_owner);
    free(face_lf);
    return face_count;
}

int mke2f(int* e2f, const int* e2n, const int* local_faces, int ne, int M, int nnf, int nfe) 
{
    int max_faces = ne * nfe;
    int* face_nodes = (int*)malloc(sizeof(int) * nnf * max_faces); // stores sorted face nodes
    int count = 0;

    for (int elem = 0; elem < ne; ++elem) {
        for (int lf = 0; lf < nfe; ++lf) {
            int temp_face[8];
            for (int i = 0; i < nnf; ++i) {
                int local_node = local_faces[i + lf * nnf];
                temp_face[i] = e2n[elem * M + local_node];
            }

            sort_face_nodes(temp_face, nnf);

            // Check if face already exists
            int match = -1;
            for (int f = 0; f < count; ++f) {
                if (equal_face_nodes(temp_face, &face_nodes[f * nnf], nnf)) {
                    match = f;
                    break;
                }
            }

            if (match >= 0) {
                e2f[elem * nfe + lf] = match;
            } else {
                // New face
                for (int i = 0; i < nnf; ++i)
                    face_nodes[count * nnf + i] = temp_face[i];

                e2f[elem * nfe + lf] = count;
                ++count;
            }
        }
    }

    free(face_nodes);
    return count;
}

void mke2e(int* e2e, const int* f2e, int nf, int ne, int nfe) {
    // Initialize e2e to -1
    for (int i = 0; i < ne * nfe; ++i)
        e2e[i] = -1;

    // Loop over all faces
    for (int i = 0; i < nf; ++i) {
        int e1 = f2e[0 + i * 4];
        int l1 = f2e[1 + i * 4];
        int e2 = f2e[2 + i * 4];
        int l2 = f2e[3 + i * 4];

        if (e1 >= 0 && e2 >= 0) {
            e2e[e1 * nfe + l1] = e2;
            e2e[e2 * nfe + l2] = e1;
        }
    }
}

void mke2f(int* e2f, const int* f2e, int nf, int nfe, int ne) 
{
    for (int i = 0; i < nfe * ne; ++i) {
        e2f[i] = -1;
    }

    for (int i = 0; i < nf; ++i) {
        int e1 = f2e[i * 4 + 0];  // f2e(:, i) �
        int l1 = f2e[i * 4 + 1];
        int e2 = f2e[i * 4 + 2];
        int l2 = f2e[i * 4 + 3];

        e2f[l1 + e1 * nfe] = i;  // e2f(l1, e1) = i

        if (e2 >= 0) {
            e2f[l2 + e2 * nfe] = i;  // e2f(l2, e2) = i
        }
    }
}

void mkf2f(int* f2f, int* f2l, const int* f2e, const int* e2f, int nf, int nfe, int ne) 
{
    int nbf = 2 * (nfe - 1);  // number of neighboring faces per face

    // Initialize f2f and f2l
    for (int i = 0; i < nbf * nf; ++i) {
        f2f[i] = -1;
        f2l[i] = -1;
    }

    for (int i = 0; i < nf; ++i) {
        int e1 = f2e[0 + 4 * i];  // element 1
        int l1 = f2e[1 + 4 * i];  // local face id on e1
        int e2 = f2e[2 + 4 * i];  // element 2 (neighbor)
        int l2 = f2e[3 + 4 * i];  // local face id on e2

        int k = 0;

        // Loop over all local faces on e1 except l1
        for (int l = 0; l < nfe; ++l) {
            if (l != l1) {
                int j = e2f[l + e1 * nfe];   // e2f(l, e1)
                f2f[k + i * nbf] = j;
                f2l[k + i * nbf] = l;
                ++k;
            }
        }

        // If interior face (e2 >= 0), repeat for second element
        if (e2 >= 0) {
            for (int l = 0; l < nfe; ++l) {
                if (l != l2) {
                    int j = e2f[l + e2 * nfe];   // e2f(l, e2)
                    f2f[k + i * nbf] = j;
                    f2l[k + i * nbf] = l;
                    ++k;
                }
            }
        }
    }
}

void mke2e(int* e2e, const int* f2e, const int* e2f, int nfe, int ne) {
    // f2e: 4 x nf, column-major, 0-based
    // e2f: nfe x ne, 0-based, index into columns of f2e (or -1 for boundary)
    // e2e: nfe x ne, output, 0-based, neighbor element index (or -1 for boundary)
    // nfe: number of faces per element
    // ne:  number of elements

    for (int i = 0; i < ne; ++i) {
        for (int j = 0; j < nfe; ++j) {
            int k = e2f[j + nfe * i]; // k: face index in f2e
            if (k < 0) {
                e2e[j + nfe * i] = -1; // boundary face
                continue;
            }
            int e1 = f2e[0 + 4 * k];
            int e2 = f2e[2 + 4 * k];
            if (e1 == i) {
                e2e[j + nfe * i] = e2; // neighbor
            } else if (e2 == i) {
                e2e[j + nfe * i] = e1; // neighbor
            } else {
                printf("Error: something wrong at element %d, local face %d\n", i, j);
            }
        }
    }
}

int mkfelem(int* felem, const int* e2f, const int* elem, int nfe, int ne) 
{
    // e2f: nfe x total_elements, 0-based, column-major
    // elem: indices of elements of interest, length ne
    // felem: output, at least nfe*ne size
    // n: output pointer, stores number of unique faces
    int n = 0;
    for (int i = 0; i < ne; ++i) {
        int el = elem[i];
        for (int j = 0; j < nfe; ++j) {
            int f = e2f[j + nfe * el]; // get face index
            int found = 0;
            for (int k = 0; k < n; ++k) {
                if (felem[k] == f) {
                    found = 1;
                    break;
                }
            }
            if (!found) {
                felem[n] = f;
                n++;
            }
        }
    }
    return n;
}

// f2e: 4 x nf, column-major, 0-based
// felem: length n, list of face indices (0-based) of interest
// elem: length ne, element indices (0-based) of interest
// n: number of faces of interest
// ne: number of elements of interest
// Output: f2eelem: 4 x n, column-major, indices remapped to [0,ne-1] if present in elem, otherwise set to 0
void mkf2eelem(int* f2eelem, const int* f2e, const int* felem, const int* elem, int n, int ne) 
{
    // Step 1: Gather relevant columns from f2e into f2eelem
    for (int i = 0; i < n; ++i) {
        for (int r = 0; r < 4; ++r)
            f2eelem[r + 4 * i] = f2e[r + 4 * felem[i]];
    }

    // Step 2: Remap element indices in f2eelem to local indices in elem (if present), otherwise handle boundaries
    for (int i = 0; i < n; ++i) {
        int e1 = f2eelem[0 + 4 * i];
        int e2 = f2eelem[2 + 4 * i];

        int found1 = 0, found2 = 0, idx1 = 0, idx2 = 0;
        for (int j = 0; j < ne; ++j) {
            if (elem[j] == e1) { found1 = 1; idx1 = j; break; }
        }
        for (int j = 0; j < ne; ++j) {
            if (elem[j] == e2) { found2 = 1; idx2 = j; break; }
        }

        if (found1) f2eelem[0 + 4 * i] = idx1;
        if (found2) f2eelem[2 + 4 * i] = idx2;

        if (!found1 && found2) {
            // Only e2 found: promote e2 info to e1 position, zero-out e2
            f2eelem[0 + 4 * i] = f2eelem[2 + 4 * i];
            f2eelem[1 + 4 * i] = f2eelem[3 + 4 * i];
            f2eelem[2 + 4 * i] = -1;
            f2eelem[3 + 4 * i] = -1;
        }
        if (!found2) {
            f2eelem[2 + 4 * i] = -1;
            f2eelem[3 + 4 * i] = -1;
        }
    }
}

void faceindex(int *in1, int *in2, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;    
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%ndf;
        int j = (idx-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int k2 = facecon[2*m+1];
        int m1 = k1%npe;
        int m2 = k2%npe;
        int n1 = (k1-m1)/npe;
        int n2 = (k2-m2)/npe;          
        in1[idx] = m1+j*npe+n1*npe*nc;
        in2[idx] = m2+j*npe+n2*npe*nc;
    }                            
}

void faceperm(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
{
    int N = 0;
    for (int j=0; j<nbf; j++) {
        int f1 = fblks[3*j]-1;
        int f2 = fblks[3*j+1];      
        int nf = f2-f1;
        int ndf = npf*nf;        
        faceindex(&ind1[N], &ind2[N], facecon, npf, ncu, npe, nc, f1, f2);
        indpts[j] = N;
        N = N + ndf*ncu;    
    }        
    indpts[nbf] = N;
}

void faceindex1(int *in1, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;    
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%ndf;
        int j = (idx-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int m1 = k1%npe;
        int n1 = (k1-m1)/npe;
        in1[idx] = m1+j*npe+n1*npe*nc;
    }                            
}

void faceperm1(int *ind1, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
{
    int N = 0;
    for (int j=0; j<nbf; j++) {
        int f1 = fblks[3*j]-1;
        int f2 = fblks[3*j+1];      
        int nf = f2-f1;
        int ndf = npf*nf;        
        faceindex1(&ind1[N], facecon, npf, ncu, npe, nc, f1, f2);
        indpts[j] = N;
        N = N + ndf*ncu;    
    }        
    indpts[nbf] = N;
}

void elemindex(int *ind, int npe, int nc, int ncu, int e1, int e2)
{        
    int nn = npe*(e2-e1);
    int N = nn*ncu;
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;   // [0, npe*ne]
        int j = (idx-i)/nn; // [0, ncu]
        int k = i%npe;  // [0, npe]
        int e = (i-k)/npe+e1;        
        ind[idx] = k+j*npe+e*npe*nc;        
    }        
}

void elemperm(int *ind, int *indpts, int *eblks, int npe, int nc, int ncu, int nbe)
{
    int N=0;
    for (int j=0; j<nbe; j++) {
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];    
        int nn = npe*(e2-e1);
        elemindex(&ind[N], npe, nc, ncu, e1, e2);
        indpts[j] = N;
        N = N + nn*ncu;
    }      
    indpts[nbe] = N;
}

int getinterfacefaces(int *bf, int *eblks, int nbe, int nfe, int ibinterface)
{ 
    int nintfaces = 0;
    for (int j=0; j<nbe; j++) { // loop over each chunk
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];               
        for (int e=e1; e<e2; e++) {  // loop over each element in a chunk           
          for (int k=0; k<nfe; k++) { // loop over each local face               
              if (bf[k + nfe*e] == ibinterface) {
                nintfaces += 1;
              }                  
          }
        }  
    }      
    
    return nintfaces;
}

void getinterfacefaces(int *intfaces, int *bf, int *eblks, int nbe, int nfe, int ibinterface, int nintfaces)
{ 
    for (int i=0; i<nintfaces; i++) intfaces[i] = 0;

    for (int j=0; j<nbe; j++) { // loop over each chunk
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];       
        int m = 0;
        for (int e=e1; e<e2; e++) {  // loop over each element in a chunk           
          for (int k=0; k<nfe; k++) { // loop over each local face               
              if (bf[k + nfe*e] == ibinterface) {
                intfaces[m] = k + nfe*e;     // local face on that element 
                m += 1;
              }                  
          }
        }  
    }                     
}


void getboundaryfaces(int *numbf, int *boufaces, int *bf, int *eblks, int nbe, int nfe, int maxbc, int nboufaces)
{
    for (int i=0; i<1+maxbc*nbe; i++) numbf[i] = 0;
    for (int i=0; i<nboufaces; i++) boufaces[i] = 0;

    for (int j=0; j<nbe; j++) { // loop over each chunk
        int e1 = eblks[3*j]-1;
        int e2 = eblks[3*j+1];            
        for (int e=e1; e<e2; e++) {  // loop over each element in a chunk           
            for (int k=0; k<nfe; k++) { // loop over each local face
                int ib = bf[k + nfe*e];                
                if (ib > 0) numbf[ib + maxbc*j] += 1; // boundary face
            }
        }          
    }                     

    // accumulative sum of numbf
    for (int j=0; j<nbe; j++) { // loop over each chunk
        for (int k=0; k<maxbc; k++) { // loop over each boundary condition
            int n = k + maxbc*j;              
            numbf[n+1] += numbf[n];
        }
    }      
          
    for (int j=0; j<nbe; j++) { // loop over each chunk
        for (int k=0; k<maxbc; k++) { // loop over each boundary condition
            int n = k + maxbc*j;
            int start = numbf[n];
            int nfaces = numbf[n+1] - start; // number of boundary faces for condition k             
            if (nfaces > 0) { // if the number of faces is not zero
              int idx = 0;
              int e1 = eblks[3*j]-1;
              int e2 = eblks[3*j+1];            
              for (int e=e1; e<e2; e++) {  // loop over each element in a chunk                    
                  for (int l=0; l<nfe; l++) { // loop over each local face
                      int ib = bf[l + nfe*e];  // boundary condition ib
                      if (ib == k+1) { // if boundary condition ib match 
                        boufaces[start + idx] = l + nfe*(e-e1); // store the boundary face (l,e)
                        idx++;
                      }
                  }
              }                                  
            }
        }
    }      
}

int getsubdomaininterfaces(const int* f2e, int ne1, int nf)
{
    int n = 0; 
    for (int j=0; j<nf; j++) { // loop over each face
        int e1 = f2e[4*j+0];   // 1st element sharing the face 
        int e2 = f2e[4*j+2];   // 2nd element sharing the face  
        if ((e1 >= ne1) || (e2 >= ne1)) { // interface faces
          n += 1;
        }
    }    
    return n;
}

int getsubdomaininterfaces(int *interface, const int* f2e, int ne1, int nf)
{
    int n = 0; 
    for (int j=0; j<nf; j++) { // loop over each face
        int e1 = f2e[4*j+0];   // 1st element sharing the face 
        int e2 = f2e[4*j+2];   // 2nd element sharing the face  
        if ((e1 >= ne1) || (e2 >= ne1)) { // interface faces
          interface[n] = j;
          n += 1;
        }
    }    
    return n;
}

void pathreorder(const int* epath, int nep, const int* e2f, int nfe,
                 int* fpath, int* lpath, int* fintf, int* lintf) {
    for (int i = 0; i < nep - 1; i++) {
        int e1 = epath[i];
        int e2 = epath[i + 1];
        const int* f1 = &e2f[e1 * nfe];
        const int* f2 = &e2f[e2 * nfe];

        int match = 0, j = -1, k = -1;
        for (int jj = 0; jj < nfe && !match; jj++) {
            for (int kk = 0; kk < nfe; kk++) {
                if (f1[jj] == f2[kk]) {
                    j = jj;
                    k = kk;
                    match = 1;
                    break;
                }
            }
        }

        if (i == 0) {
            int m = 0;
            if (nfe == 8) m = (j % 2 == 0) ? j - 1 : j + 1;
            else if (nfe == 4) m = (j == 0) ? 2 : (j == 2) ? 0 : (j == 1) ? 3 : 1;
            else if (nfe == 3) m = (j == 0) ? 2 : (j == 1) ? 0 : 1;

            fpath[i * 2 + 0] = f1[m];
            fpath[i * 2 + 1] = f1[j];
            fpath[(i + 1) * 2 + 0] = f2[k];
            lpath[i * 2 + 0] = m;
            lpath[i * 2 + 1] = j;
            lpath[(i + 1) * 2 + 0] = k;
        }

        if (i == nep - 2) {
            int m = 0;
            if (nfe == 8) m = (k % 2 == 0) ? k - 1 : k + 1;
            else if (nfe == 4) m = (k == 0) ? 2 : (k == 2) ? 0 : (k == 1) ? 3 : 1;
            else if (nfe == 3) m = (k == 0) ? 2 : (k == 1) ? 0 : 1;

            fpath[(nep - 1) * 2 + 0] = f2[k];
            fpath[(nep - 1) * 2 + 1] = f2[m];
            fpath[(nep - 2) * 2 + 1] = f1[j];
            lpath[(nep - 1) * 2 + 0] = k;
            lpath[(nep - 1) * 2 + 1] = m;
            lpath[(nep - 2) * 2 + 1] = j;
        }

        if (i != 0 && i != nep - 2) {
            fpath[i * 2 + 1] = f1[j];
            fpath[(i + 1) * 2 + 0] = f2[k];
            lpath[i * 2 + 1] = j;
            lpath[(i + 1) * 2 + 0] = k;
        }
    }

    for (int i = 0; i < nep; i++) {
        int j = lpath[i * 2 + 0];
        int m = lpath[i * 2 + 1];
        const int* fi = &e2f[epath[i] * nfe];
        int n = 0;
        for (int l = 0; l < nfe; l++) {
            if (l != m && l != j) {
                fintf[i * (nfe - 2) + n] = fi[l];
                lintf[i * (nfe - 2) + n] = l;
                n++;
            }
        }
    }
}

void pathreordering(int* fpath, int* lpath, int* fintf, int* lintf, const int* epath, const int* e2f, int nfe, int npaths, int nep) 
{
    int* tmp_epath = (int*)malloc(nep * sizeof(int));
    int* tmp_fpath = (int*)malloc(2 * nep * sizeof(int));
    int* tmp_lpath = (int*)malloc(2 * nep * sizeof(int));
    int* tmp_fintf = (int*)malloc((nfe - 2) * nep * sizeof(int));
    int* tmp_lintf = (int*)malloc((nfe - 2) * nep * sizeof(int));

    for (int i = 0; i < npaths; i++) {
        for (int j = 0; j < nep; j++) tmp_epath[j] = epath[i + j * npaths];
        pathreorder(tmp_epath, nep, e2f, nfe,
                    tmp_fpath, tmp_lpath, tmp_fintf, tmp_lintf);

        for (int j = 0; j < nep; j++) {
            fpath[(i + j * npaths) * 2 + 0] = tmp_fpath[j * 2 + 0];
            fpath[(i + j * npaths) * 2 + 1] = tmp_fpath[j * 2 + 1];
            lpath[(i + j * npaths) * 2 + 0] = tmp_lpath[j * 2 + 0];
            lpath[(i + j * npaths) * 2 + 1] = tmp_lpath[j * 2 + 1];
            for (int k = 0; k < nfe - 2; k++) {
                fintf[(i + j * npaths) * (nfe - 2) + k] = tmp_fintf[j * (nfe - 2) + k];
                lintf[(i + j * npaths) * (nfe - 2) + k] = tmp_lintf[j * (nfe - 2) + k];
            }
        }
    }

    free(tmp_epath);
    free(tmp_fpath);
    free(tmp_lpath);
    free(tmp_fintf);
    free(tmp_lintf);
}


void simple_bubble_sort(int* b, int* ind, const int* a, int n) 
{
    // Initialize b and ind
    for (int i = 0; i < n; ++i) {
        b[i] = a[i];
        ind[i] = i;
    }
    // Bubble sort
    for (int i = 0; i < n-1; ++i) {
        for (int j = 0; j < n-i-1; ++j) {
            if (b[j] > b[j+1]) {
                // Swap values
                int tmp = b[j];
                b[j] = b[j+1];
                b[j+1] = tmp;
                // Swap indices
                int tmp_idx = ind[j];
                ind[j] = ind[j+1];
                ind[j+1] = tmp_idx;
            }
        }
    }
}

void printintarray(Int* a, Int m, Int n)
{
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void sortinteriorfaces(int* f2e, int nfe, int ne, int nf, int nf0) 
{
    int nfsorted = nf0;
    int* e2f = (int*)malloc(nfe * ne * sizeof(int));
    int nbf = 2 * (nfe - 1);

    int* f2f = (int*)malloc(nbf * nf * sizeof(int));
    int* f2l = (int*)malloc(nbf * nf * sizeof(int));
    int* unsorted = (int*)malloc(nf * sizeof(int));         // largest possible
    int* unsortedcount = (int*)malloc(nf * sizeof(int));
    int* a = (int*)malloc(nf * sizeof(int));                // for permutation indices
    int* unsortedcount_sorted = (int*)malloc(nf * sizeof(int));

    while (nfsorted < nf) {
        // Build e2f and f2f
        mke2f(e2f, f2e, nf, nfe, ne);           // e2f: nfe x ne
        mkf2f(f2f, f2l, f2e, e2f, nf, nfe, ne); // f2f: nbf x nf

        int nfunsorted = nf - nfsorted;
        // Build list of unsorted face indices
        for (int i = 0; i < nfunsorted; ++i)
            unsorted[i] = nfsorted + i;
        
        // Compute unsortedcount
        for (int i = 0; i < nfunsorted; ++i) {
            int fi = unsorted[i];
            unsortedcount[i] = 0;
            for (int j = 0; j < nbf; ++j) {
                int fnei = f2f[j + nbf * fi];
                if (fnei >= nfsorted) {
                    unsortedcount[i]++;
                }
            }
        }        
        
        //printintarray(unsortedcount, 1, nfunsorted);
                
        // Sort unsortedcount (in ascending order), get permutation a
        simple_bubble_sort(unsortedcount_sorted, a, unsortedcount, nfunsorted);
        //void simple_bubble_sort(int* b, int* ind, const int* a, int n) 
        
        // Reorder columns of f2e for unsorted faces
        // Save original columns for these faces
        int* tmpcols = (int*)malloc(4 * nfunsorted * sizeof(int));
        for (int i = 0; i < nfunsorted; ++i) {
            for (int r = 0; r < 4; ++r) {
                tmpcols[r + 4 * i] = f2e[r + 4 * unsorted[a[i]]];
            }
        }
        // Now assign back to f2e in sorted order
        for (int i = 0; i < nfunsorted; ++i) {
            int newidx = unsorted[i];
            for (int r = 0; r < 4; ++r) {
                f2e[r + 4 * newidx] = tmpcols[r + 4 * i];
            }
        }
        free(tmpcols);
        
//        printintarray(f2e, 4, nf);

        // Update nfsorted
        int inc = 0;
        for (int i = 0; i < nfunsorted; ++i)
            if (unsortedcount_sorted[i] < nbf)
                inc++;
        nfsorted += inc;        
    }

    free(e2f);
    free(f2f);
    free(f2l);
    free(unsorted);
    free(unsortedcount);
    free(a);
    free(unsortedcount_sorted);
}

// Sorts array by recurrence, returns b (sorted array), e (indices mapping b to a)
void sortrecurrence(int* b, int* e, const int* a, int n) 
{
    int* recurrence = (int*)malloc(n * sizeof(int));

    // Step 1: Count occurrences for each entry
    for (int i = 0; i < n; ++i) {
        int count = 0;
        for (int j = 0; j < n; ++j) {
            if (a[i] == a[j])
                count++;
        }
        recurrence[i] = count;
    }

    // Find maxrec
    int maxrec = 0;
    for (int i = 0; i < n; ++i) {
        if (recurrence[i] > maxrec)
            maxrec = recurrence[i];
    }

    int* c = (int*)malloc(n * sizeof(int));
    int* d = (int*)malloc(n * sizeof(int));
    int* c_sorted = (int*)malloc(n * sizeof(int));
    int* ind = (int*)malloc(n * sizeof(int));
    int m = 0;
    for (int k = 1; k <= maxrec; ++k) {
        int j = 0;
        for (int i = 0; i < n; ++i) {
            if (recurrence[i] == k) {
                c[j] = a[i];
                d[j] = i;
                j++;
            }
        }
        if (j > 0) {
            // Sort c and reorder d accordingly
            simple_bubble_sort(c_sorted, ind, c, j);
            //void simple_bubble_sort(int* b, int* ind, const int* a, int n) 
            
            // Apply permutation to d
            int* d_sorted = (int*)malloc(j * sizeof(int));
            for (int ii = 0; ii < j; ++ii)
                d_sorted[ii] = d[ind[ii]];
            // Copy to output
            for (int ii = 0; ii < j; ++ii) {
                b[m] = c_sorted[ii];
                e[m] = d_sorted[ii];
                m++;
            }
            free(d_sorted);
        }
        if (m >= n)
            break;
    }

    // Check that b equals a[e]
    for (int i = 0; i < n; ++i) {
        if (b[i] != a[e[i]]) {
            fprintf(stderr, "something wrong at i = %d: b[i]=%d, a[e[i]]=%d\n", i, b[i], a[e[i]]);
            exit(EXIT_FAILURE);
        }
    }

    free(recurrence);
    free(c);
    free(d);
    free(c_sorted);
    free(ind);
}

// Helper to find indices where the last row of f2e is -1 (boundary faces)
int find_boundary_faces(int* ind0, const int* f2e, int nf) {
    int count = 0;
    for (int i = 0; i < nf; ++i) {
        if (f2e[3 + 4*i] == -1) { // Last row, 0-based
            ind0[count++] = i;
        }
    }
    return count; // number of boundary faces found
}

// Helper to compute setdiff between 0..nf-1 and ind0[0..nind0-1]
int setdiff(int* ind1, const int* ind0, int nind0, int nf) {
    int count = 0;
    for (int i = 0; i < nf; ++i) {
        int found = 0;
        for (int j = 0; j < nind0; ++j)
            if (i == ind0[j]) { found = 1; break; }
        if (!found) ind1[count++] = i;
    }
    return count; // number of interior faces
}

void facereorder(int* f2e, int nf) 
{
    // Step 1: find boundary face indices (last row == 0)
    int* ind0 = (int*)malloc(nf * sizeof(int));
    int nind0 = find_boundary_faces(ind0, f2e, nf);
    
    // Step 2: sort boundary faces by recurrence of first row
    int* firstrow = (int*)malloc(nind0 * sizeof(int));
    for (int i = 0; i < nind0; ++i)
        firstrow[i] = f2e[0 + 4*ind0[i]];
    int* tmpb = (int*)malloc(nind0 * sizeof(int));
    int* ii = (int*)malloc(nind0 * sizeof(int));
    sortrecurrence(tmpb, ii, firstrow, nind0); // ii: permutation
    // Reorder ind0 according to ii
    int* ind0_sorted = (int*)malloc(nind0 * sizeof(int));
    for (int i = 0; i < nind0; ++i)
        ind0_sorted[i] = ind0[ii[i]];
    
    for (int i = 0; i < nind0; ++i)
      ind0[i] = ind0_sorted[i];
    
    // Step 3: get indices of interior faces
    int* ind1 = (int*)malloc(nf * sizeof(int));
    int nind1 = setdiff(ind1, ind0, nind0, nf);

    // Step 4: build the reordered face list
    int* order = (int*)malloc(nf * sizeof(int));
    for (int i = 0; i < nind0; ++i)
        order[i] = ind0[i];
    for (int i = 0; i < nind1; ++i)
        order[nind0 + i] = ind1[i];
    
    //printintarray(order, 1, nf);

    // Step 5: permute columns of f2e
    int* f2e_new = (int*)malloc(4 * nf * sizeof(int));
    for (int i = 0; i < nf; ++i)
        for (int j = 0; j < 4; ++j)
            f2e_new[j + 4*i] = f2e[j + 4*order[i]];
    memcpy(f2e, f2e_new, 4*nf*sizeof(int));
    free(f2e_new);
    
    //printintarray(f2e, 4, nf);

    // Step 6: find ne and nfe
    int ne1 = 0, ne2 = 0, nfe1 = 0, nfe2 = 0;
    for (int i = 0; i < nf; ++i) {
        if (f2e[0 + 4*i] > ne1) ne1 = f2e[0 + 4*i];
        if (f2e[2 + 4*i] > ne2) ne2 = f2e[2 + 4*i];
        if (f2e[1 + 4*i] > nfe1) nfe1 = f2e[1 + 4*i];
        if (f2e[3 + 4*i] > nfe2) nfe2 = f2e[3 + 4*i];
    }
    int ne = (ne1 > ne2) ? ne1 : ne2;
    int nfe = (nfe1 > nfe2) ? nfe1 : nfe2;
    ne += 1;
    nfe += 1;
    
    // Step 7: sort interior faces
    sortinteriorfaces(f2e, nfe, ne, nf, nind0);

    //printintarray(f2e, 4, nf);
    
    // Free temporary arrays
    free(ind0);
    free(ind0_sorted);
    free(ind1);
    free(order);
    free(firstrow);
    free(tmpb);
    free(ii);
}

void bubble_sort(int* a, int n) {
    for (int i = 0; i < n-1; ++i)
        for (int j = 0; j < n-i-1; ++j)
            if (a[j] > a[j+1]) {
                int tmp = a[j]; a[j] = a[j+1]; a[j+1] = tmp;
            }
}

void crs_array(int** row_ptr_out, int** col_ind_out, const int* f2e, int nfe, int ne, int nf) 
{
    // Step 1: Compute e2f and f2f
    int* e2f = (int*)malloc(nfe * ne * sizeof(int));
    mke2f(e2f, f2e, nf, nfe, ne);

    int nbf = 2 * (nfe - 1);
    int* f2f = (int*)malloc(nbf * nf * sizeof(int));
    int* f2l = (int*)malloc(nbf * nf * sizeof(int));
    mkf2f(f2f, f2l, f2e, e2f, nf, nfe, ne);

    // Step 2: row_ptr
    int* row_ptr = (int*)malloc((nf + 1) * sizeof(int));
    row_ptr[0] = 0;
    for (int i = 0; i < nf; ++i) {
        int n = 1;
        for (int j = 0; j < nbf; ++j)
            if (f2f[j + nbf * i] >= 0)
                n++;
        row_ptr[i + 1] = row_ptr[i] + n;
    }

    // Step 3: col_ind
    int nnz = row_ptr[nf];
    int* col_ind = (int*)malloc(nnz * sizeof(int));
    for (int i = 0; i < nf; ++i) {
        // Count how many entries f2f(:,i) > 0
        int cnt = 0;
        for (int j = 0; j < nbf; ++j)
            if (f2f[j + nbf * i] >= 0)
                cnt++;
        // Gather these values
        int* neighbors = (int*)malloc(cnt * sizeof(int));
        int k = 0;
        for (int j = 0; j < nbf; ++j)
            if (f2f[j + nbf * i] >= 0)
                neighbors[k++] = f2f[j + nbf * i];
        // Sort the neighbors
        bubble_sort(neighbors, cnt);
        // Assign col_ind: first the diagonal (i), then sorted neighbors
        int base = row_ptr[i];
        col_ind[base] = i;
        for (int j = 0; j < cnt; ++j)
            col_ind[base + 1 + j] = neighbors[j];
        free(neighbors);
    }

    // Return results
    *row_ptr_out = row_ptr;
    *col_ind_out = col_ind;

    free(e2f);
    free(f2f);
    free(f2l);
}

void mkface(int* face, const int* e2f, const int* elem, const int* f2eelem, int nb, int nf, int nfe) 
{
    for (int n = 0; n < nb; ++n) {
        for (int i = 0; i < nf; ++i) {
            int e = f2eelem[0 + 4*i]; // 0-based
            int l = f2eelem[1 + 4*i];
            int local_e = elem[n + nb*e]; // local node index for this element
            face[n + nb*i] = e2f[l + nfe*local_e];
        }
    }
}

int crs_faceordering(
    int** row_ptr_out,
    int** col_ind_out,
    int** face_out,
    int** f2eelem_out,    
    const int* elem,  // nb x ne, column-major    
    const int* f2e,   // 4 x nf, column-major
    int nb, int ne, int nfe, int nf            
) {
    // Step 1: find mesh parameters (0-based!)  
    int ne1=0, ne2=0;
    for (int i = 0; i < nf; ++i) if (f2e[0 + 4*i] > ne1) ne1 = f2e[0 + 4*i];    
    for (int i = 0; i < nf; ++i) if (f2e[2 + 4*i] > ne2) ne2 = f2e[0 + 4*i];    
    int ne_tot = (ne1 > ne2) ? ne1 : ne2;
    ne_tot = ne_tot+1;
    
    // Step 2: build e2f
    int* e2f = (int*)malloc(nfe * ne_tot * sizeof(int));
    mke2f(e2f, f2e, nf, nfe, ne_tot);
    
    // Step 3: mkfelem
    int* elem_row0 = (int*)malloc(ne * sizeof(int));
    for (int i = 0; i < ne; ++i) elem_row0[i] = elem[0 + nb * i];
    int* felem = (int*)malloc(nfe * ne * sizeof(int)); // upper bound
    int nfelem = mkfelem(felem, e2f, elem_row0, nfe, ne);
    //int mkfelem(int* felem, const int* e2f, const int* elem, int nfe, int ne) 
    
    //printintarray(elem_row0, 1, ne);
    //printintarray(felem, 1, nfelem);
        
    // Step 4: mkf2eelem
    int* f2eelem = (int*)malloc(4 * nfelem * sizeof(int));
    mkf2eelem(f2eelem, f2e, felem, elem_row0, nfelem, ne);
    // void mkf2eelem(int* f2eelem, const int* f2e, const int* felem, const int* elem, int n, int ne) 
    
    //printintarray(f2eelem, 4, nfelem);
        
    // Step 5: facereorder
    facereorder(f2eelem, nfelem);
    // void facereorder(int* f2e, int nf) 
        
    // Step 6: crs_array
    int* row_ptr = NULL;
    int* col_ind = NULL;
    crs_array(&row_ptr, &col_ind, f2eelem, nfe, ne, nfelem);
    // void crs_array(int** row_ptr_out, int** col_ind_out, const int* f2e, int nfe, int ne, int nf) 
    
    // Step 7: mkface
    int* face = (int*)malloc(nb * nfelem * sizeof(int));
    mkface(face, e2f, elem, f2eelem, nb, nfelem, nfe);
    // void mkface(int* face, const int* e2f, const int* elem, const int* f2eelem, int nb, int nf, int nfe) 
        
    // Output
    *row_ptr_out = row_ptr;
    *col_ind_out = col_ind;
    *face_out = face;
    *f2eelem_out = f2eelem;
    
    // Free intermediate
    free(e2f);
    free(felem);
    free(elem_row0);
        
    return nfelem;
}

int gridpartition2d(int **elemblocks_out, int Nx, int Ny, int Mx, int My, int opt)
{
    // Check input
    if (Nx % Mx != 0 || Ny % My != 0)
    {
        printf("Error: Nx and Ny must be divisible by Mx and My\n");
        *elemblocks_out = NULL;
        return 0;
    }
    int nblocks_x = Nx / Mx;
    int nblocks_y = Ny / My;
    int nblocks = nblocks_x * nblocks_y;
    int blocksize = Mx * My;

    // Create grid
    int *grid = (int *)malloc(Nx * Ny * sizeof(int));
    int *elemblocks = (int *)malloc(nblocks * blocksize * sizeof(int));

    if (opt == 0)
    {
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i)
                grid[i + Nx * j] = i + Nx * j;

        int blockidx = 0;
        for (int j = 0; j < nblocks_y; ++j)
        {
            for (int i = 0; i < nblocks_x; ++i)
            {
                int cnt = 0;
                int x_start = i * Mx;
                int y_start = j * My;
                for (int jj = 0; jj < My; ++jj)
                {
                    for (int ii = 0; ii < Mx; ++ii)
                    {
                        int gx = x_start + ii;
                        int gy = y_start + jj;
                        elemblocks[blockidx * blocksize + cnt++] = grid[gx + Nx * gy];
                    }
                }
                blockidx++;
            }
        }
    }
    else
    {
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                grid[i + Nx * j] = j + Ny * i;

        int blockidx = 0;
        for (int i = 0; i < nblocks_x; ++i)
        {
            for (int j = 0; j < nblocks_y; ++j)
            {
                int cnt = 0;
                int x_start = i * Mx;
                int y_start = j * My;
                for (int jj = 0; jj < My; ++jj)
                {
                    for (int ii = 0; ii < Mx; ++ii)
                    {
                        int gx = x_start + ii;
                        int gy = y_start + jj;
                        elemblocks[blockidx * blocksize + cnt++] = grid[gx + Nx * gy];
                    }
                }
                blockidx++;
            }
        }
    }

    int *elemblocks_t = (int *)malloc(nblocks * blocksize * sizeof(int));
    for (int i=0; i<blocksize; i++)
      for (int j=0; j<nblocks; j++)
        elemblocks_t[j + nblocks*i] = elemblocks[i + blocksize*j];
    
    free(grid);
    free(elemblocks);
    *elemblocks_out = elemblocks_t;
    return nblocks;    
}

int gridpartition3d(int** elemblocks_out, int Nx, int Ny, int Nz, int Mx, int My, int Mz, int opt) 
{
    // Check divisibility
    if (Nx % Mx != 0 || Ny % My != 0 || Nz % Mz != 0) {
        printf("Error: Grid dimensions must be divisible by block sizes.\n");
        *elemblocks_out = NULL;
        return 0;
    }

    int nblocks_x = Nx / Mx;
    int nblocks_y = Ny / My;
    int nblocks_z = Nz / Mz;
    int nblocks = nblocks_x * nblocks_y * nblocks_z;
    int blocksize = Mx * My * Mz;

    // Allocate memory for grid and blocks (all 0-based)
    int* grid = (int*)malloc(Nx * Ny * Nz * sizeof(int));
    int* elemblocks = (int*)malloc(nblocks * blocksize * sizeof(int));

    if (opt == 0) {
        // Fill grid with lexicographical indices (i + Nx*(j) + Nx*Ny*(k))
        for (int k = 0; k < Nz; ++k)
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    grid[i + Nx * (j + Ny * k)] = i + Nx * j + Nx * Ny * k;

        int blockidx = 0;
        for (int kk = 0; kk < nblocks_z; ++kk)
            for (int jj = 0; jj < nblocks_y; ++jj)
                for (int ii = 0; ii < nblocks_x; ++ii) {
                    int cnt = 0;
                    int x_start = ii * Mx;
                    int y_start = jj * My;
                    int z_start = kk * Mz;
                    for (int k = 0; k < Mz; ++k)
                        for (int j = 0; j < My; ++j)
                            for (int i = 0; i < Mx; ++i) {
                                int gx = x_start + i;
                                int gy = y_start + j;
                                int gz = z_start + k;
                                elemblocks[blockidx * blocksize + cnt++] = grid[gx + Nx * (gy + Ny * gz)];
                            }
                    blockidx++;
                }
    } else {
        // Fill grid with "column-major" indices (k + Nz*j + Nz*Ny*i)
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j)
                for (int k = 0; k < Nz; ++k)
                    grid[i + Nx * (j + Ny * k)] = k + Nz * j + Nz * Ny * i;

        int blockidx = 0;
        for (int ii = 0; ii < nblocks_x; ++ii)
            for (int jj = 0; jj < nblocks_y; ++jj)
                for (int kk = 0; kk < nblocks_z; ++kk) {
                    int cnt = 0;
                    int x_start = ii * Mx;
                    int y_start = jj * My;
                    int z_start = kk * Mz;
                    for (int k = 0; k < Mz; ++k)
                        for (int j = 0; j < My; ++j)
                            for (int i = 0; i < Mx; ++i) {
                                int gx = x_start + i;
                                int gy = y_start + j;
                                int gz = z_start + k;
                                elemblocks[blockidx * blocksize + cnt++] = grid[gx + Nx * (gy + Ny * gz)];
                            }
                    blockidx++;
                }
    }

    int *elemblocks_t = (int *)malloc(nblocks * blocksize * sizeof(int));
    for (int i=0; i<blocksize; i++)
      for (int j=0; j<nblocks; j++)
        elemblocks_t[j + nblocks*i] = elemblocks[i + blocksize*j];
    
    free(grid);
    free(elemblocks);
    *elemblocks_out = elemblocks_t;
    return nblocks;    
}

void crs_indexingilu0(
    // Output arrays (all must be pre-allocated by the caller)
    int* ind_ii,        // [N]
    int* ind_ji,        // [M*N], M = 2*(nfe-1)
    int* ind_jl,        // [M*M*N]
    int* ind_il,        // [M*M*N]
    int* num_ji,        // [N]
    int* num_jl,        // [M*N]
    int* Lind_ji,       // [M*2*N]
    int* Uind_ji,       // [M*2*N]
    int* Lnum_ji,       // [2*N]
    int* Unum_ji,       // [3*N]        
    const int* row_ptr, // size N+1
    const int* col_ind, // CRS col indices, size nnz
    int nfe,            // num. faces per element
    int N               // number of rows/nodes
)
{
    int M = 2 * (nfe - 1);

    // 1. Main loop: ind_ii, ind_ji, ind_jl, ind_il, num_ji, num_jl
    for (int i = 0; i < N; ++i) {
        int r1 = row_ptr[i];
        int r2 = row_ptr[i+1] - 1; // inclusive
        int diag_idx = -1;
        for (int p = r1; p <= r2; ++p) {
            if (col_ind[p] == i) {
                diag_idx = p;
                break;
            }
        }
        if (diag_idx == -1) {
            printf("ILU0: missing diagonal block at row %d\n", i);
            exit(1);
        }
        ind_ii[i] = diag_idx;

        int k = 0;
        for (int p = r1; p <= r2; ++p) {
            int j = col_ind[p];
            if (j <= i) continue;

            int rj1 = row_ptr[j];
            int rj2 = row_ptr[j+1] - 1;
            int idx_ji = -1;
            for (int q = rj1; q <= rj2; ++q) {
                if (col_ind[q] == i) {
                    idx_ji = q;
                    break;
                }
            }
            if (idx_ji == -1) continue;

            // ind_ji: [k + i*M]
            ind_ji[k + i*M] = idx_ji;
            num_ji[i] = k+1;

            int m = 0;
            for (int pp = r1; pp <= r2; ++pp) {
                int ell = col_ind[pp];
                if (ell <= i) continue;

                int idx_jl = -1;
                for (int qq = rj1; qq <= rj2; ++qq) {
                    if (col_ind[qq] == ell) {
                        idx_jl = qq;
                        break;
                    }
                }
                if (idx_jl == -1) continue;

                // ind_jl, ind_il: [m + k*M + i*M*M]
                ind_jl[m + k*M + i*M*M] = idx_jl;
                ind_il[m + k*M + i*M*M] = pp;
                num_jl[k + i*M] = m+1;
                m++;
            }
            k++;
        }
    }

    // 2. Lind_ji, Lnum_ji
    for (int i = 0; i < N; ++i) {
        int rstart = row_ptr[i];
        int rend   = row_ptr[i+1] - 1;
        int k = 0;
        for (int ptr = rstart; ptr <= rend; ++ptr) {
            int j = col_ind[ptr];
            if (j < i) {
                Lind_ji[k + 0*M + i*2*M] = ptr; // (k,1,i)
                Lind_ji[k + 1*M + i*2*M] = j;   // (k,2,i)
                k++;
            }
        }
        Lnum_ji[0 + 2*i] = k;
        Lnum_ji[1 + 2*i] = 1;
        for (int l = 1; l < k; ++l) {
            if (Lind_ji[l + 0*M + i*2*M] - Lind_ji[l-1 + 0*M + i*2*M] != 1) {
                Lnum_ji[1 + 2*i] = 0;
                break;
            }
        }
    }

    // 3. Uind_ji, Unum_ji
    for (int i = N-1; i >= 0; --i) {
        int rstart = row_ptr[i];
        int rend   = row_ptr[i+1] - 1;
        int k = 0;
        for (int ptr = rstart; ptr <= rend; ++ptr) {
            int ell = col_ind[ptr];
            if (ell > i) {
                Uind_ji[k + 0*M + i*2*M] = ptr; // (k,1,i)
                Uind_ji[k + 1*M + i*2*M] = ell;// (k,2,i)
                k++;
            }
        }
        Unum_ji[0 + 3*i] = k;
        Unum_ji[1 + 3*i] = rstart;
        Unum_ji[2 + 3*i] = 1;
        for (int l = 1; l < k; ++l) {
            if (Uind_ji[l + 0*M + i*2*M] - Uind_ji[l-1 + 0*M + i*2*M] != 1) {
                Unum_ji[2 + 3*i] = 0;
                break;
            }
        }
    }
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

void gaussnodes(vector<dstype>& xgauss, vector<dstype>& wgauss,
                int pgauss, int dim, int elemtype, const std::string filename) 
{
    ifstream file(filename, ios::binary);
    if (!file) error("Error opening file: " + filename);

    // Read the file into a buffer
    file.seekg(0, ios::end);
    size_t num_bytes = file.tellg();
    file.seekg(0, ios::beg);
    size_t num_doubles = num_bytes / sizeof(dstype);
    vector<dstype> tmp(num_doubles);
    file.read(reinterpret_cast<char*>(tmp.data()), num_bytes);
    file.close();

    // Read header
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

    // Compute cumulative lengths
    vector<int> lz(total_blocks + 1, 0);
    partial_sum(sz.begin(), sz.end(), lz.begin() + 1);

    int data_start = offset + 2 * total_blocks;

    auto extract_block = [&](int i, vector<dstype>& out) {
        // Corrected zero-based indexing
        int e = elemtype + 1;
        int idx = index4D(i, e - 1, pgauss - 1, dim - 1, narrays);
        int start = lz[idx];
        int count = sz[idx];
        out.resize(count);
        copy(tmp.begin() + data_start + start,
             tmp.begin() + data_start + start + count,
             out.begin());
    };

    extract_block(0, xgauss);
    extract_block(1, wgauss);
}

#endif  

