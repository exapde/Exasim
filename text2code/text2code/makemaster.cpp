/*
 * makemaster.cpp
 * 
 * This file contains functions and structures for initializing and managing the "Master" element in high-order finite element methods,
 * including shape function generation, nodal and quadrature point management, and mesh construction utilities.
 * 
 * Main Components:
 * 
 * - struct Master:
 *     Stores all master element data, including nodal and quadrature points, shape functions, permutation indices, and element metadata.
 * 
 * - int index4D(...):
 *     Computes a column-major index for 4D arrays, used for block extraction from binary files.
 * 
 * - void masternodes(...):
 *     Reads nodal and face data for master elements from a binary file, parses header, and extracts blocks for element/face nodes and permutations.
 * 
 * - void gaussnodes(...):
 *     Reads Gauss quadrature nodes and weights from a binary file for integration over elements/faces.
 * 
 * - void pascalindex2d(...), void pascalindex3d(...):
 *     Generates Pascal triangle/tetrahedron indices for polynomial basis construction in 2D/3D.
 * 
 * - void legendre(...), void jacobi(...):
 *     Computes coefficients for Legendre and Jacobi polynomials, used in orthogonal basis construction.
 * 
 * - void polyder(...), void polyval(...), void conv(...):
 *     Polynomial utilities for differentiation, evaluation, and convolution.
 * 
 * - void koornwinder1d/2d/3d(...), void koornwinder(...):
 *     Constructs Koornwinder orthogonal polynomial basis functions and their derivatives for 1D, 2D, and 3D elements.
 * 
 * - void tensorproduct(...):
 *     Constructs tensor-product basis functions for quadrilateral/hexahedral elements.
 * 
 * - void mkshape(...):
 *     Generates nodal shape functions and their derivatives by inverting the Vandermonde matrix of orthogonal basis functions.
 * 
 * - void localbasis(...):
 *     Computes local basis functions for vertices and faces of standard elements (line, triangle, quadrilateral, tetrahedron, hexahedron).
 * 
 * - int permindex(...):
 *     Computes permutation indices for face nodes under symmetry operations for different element types and dimensions.
 * 
 * - Master initializeMaster(...):
 *     Initializes the Master structure by reading nodal/quadrature data, generating shape functions, and setting up element metadata.
 * 
 * - void writemaster(...):
 *     Writes all master element data to a binary file for later use.
 * 
 * - void buildMesh(...):
 *     Constructs mesh connectivity, boundary/periodic faces, and projects DG nodes onto curved boundaries using master element data.
 * 
 * Usage:
 *   - Used in high-order DG/finite element codes for initializing reference element data, shape functions, and mesh connectivity.
 *   - Supports both simplex (triangle/tetrahedron) and tensor-product (quadrilateral/hexahedron) elements in 1D, 2D, and 3D.
 * 
 * Dependencies:
 *   - Requires LAPACK routines (DGETRF, DGETRI, DGEMM) for matrix inversion and multiplication.
 *   - Relies on external binary files for nodal and quadrature data.
 * 
 * Author: Exasim Project
 * License: MIT
 */

#ifndef __MAKEMASTER
#define __MAKEMASTER

struct Master 
{    
    vector<double> xpe, gpe, gwe, xpf, gpf, gwf; // shapegwdotshapeg, shapfgwdotshapfg;
    vector<double> shapeg, shapegt, shapfg, shapfgt, shapent, shapfnt, shapegw, shapfgw, shapen, shapfn;
    vector<double> xp1d, gp1d, gw1d, shap1dg, shap1dgt, shap1dn, shap1dnt, shap1dgw, phielem, phiface;
    vector<int>  telem, tface, perm, permind; 
    int nd, npe, npf, nge, ngf, porder, pgauss, nfe, elemtype, nodetype, nve, nvf, np1d, ng1d, npermind;
};

int index4D(int i, int j, int k, int l, const vector<int>& shape) {
    // Column-major indexing: idx = i + j*n1 + k*n1*n2 + l*n1*n2*n3
    return i + shape[0] * (j + shape[1] * (k + shape[2] * l));
}

void masternodes(vector<double>& pelem, vector<int>& telem,
                 vector<double>& pface, vector<int>& tface,
                 vector<int>& perm, int porder, int dim, int elemtype, const std::string filename) 
{
    
    ifstream file(filename, ios::binary);
    
    if (!file) error("Error opening file: " + filename);

    // Read the full file into a vector
    file.seekg(0, ios::end);
    size_t num_bytes = file.tellg();
    file.seekg(0, ios::beg);
    size_t num_doubles = num_bytes / sizeof(double);

    vector<double> tmp(num_doubles);
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

//     printf("ndims = %d\n", ndims);
//     printf("total_blocks = %d\n", total_blocks);
//     printf("offset = %d\n", offset);
//     printf("data_start = %d\n", data_start);
//     print2iarray(narrays.data(), 1, ndims);
//     print2iarray(sz1.data(), 1, total_blocks);
//     print2iarray(sz2.data(), 1, total_blocks);
//     print2iarray(sz.data(), 1, total_blocks);
//     print2iarray(lz.data(), 1, total_blocks+1);
    
    auto extract_block = [&](int i, vector<double>& out) {
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

    vector<double>  telemd, tfaced, permd;
    
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

void gaussnodes(vector<double>& xgauss, vector<double>& wgauss,
                int pgauss, int dim, int elemtype, const std::string filename) 
{
    ifstream file(filename, ios::binary);
    if (!file) error("Error opening file: " + filename);

    // Read the file into a buffer
    file.seekg(0, ios::end);
    size_t num_bytes = file.tellg();
    file.seekg(0, ios::beg);
    size_t num_doubles = num_bytes / sizeof(double);
    vector<double> tmp(num_doubles);
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

    auto extract_block = [&](int i, vector<double>& out) {
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

void pascalindex2d(int *pq, int numPoly)
{
    int i, j, counter = 0, done = 0;
    
    for (i = 0; i <= numPoly; i++) {
        for (j = 0; j <= i; j++) {
            pq[0*numPoly+counter] = i-j;
            pq[1*numPoly+counter] = j;
            counter++;
            if (counter >= numPoly) {
                done = 1;
                break;
            }
        }
        if (done == 1)
            break;
    }
}

void pascalindex3d(int *pq, int numPoly)
{
    int i, j, k, counter = 0, done = 0;
    
    for (i = 0; i <= numPoly; i++) {
        for (j = 0; j <= i; j++) {
            for (k = 0; k <= j; k++) {
                pq[0*numPoly+counter] = i-j;
                pq[1*numPoly+counter] = j-k;
                pq[2*numPoly+counter] = k;
                counter++;
                if (counter >= numPoly) {
                    done = 1;
                    break;
                }
            }
            if (done == 1)
                break;
        }
        if (done == 1)
            break;
    }
}

void legendre(double *coeff, int porder)
{
    int i;
    
    if (porder == 0)        // P_0 = 1.0
        coeff[0] = 1.0;
    else if (porder == 1) { // P_1 = x
        coeff[0] = 1.0;
        coeff[1] = 0.0;
    }
    else {                  // (n+1) * P_{n+1} = (2*n+1)*x*P_n - n*P_{n-1}
        double *coeff1 = new double[porder];
        double *coeff2 = new double[porder-1];
        
        legendre(&coeff1[0], porder-1);
        legendre(&coeff2[0], porder-2);
        
        for (i = 0; i <= porder; i++)
            coeff[i] = 0.0;
        for (i = 0; i <= porder-1; i++)
            coeff[i] += (2.0 * ((double) porder) - 1.0) * coeff1[i] / (double) porder;
        for (i = 0; i <= porder-2; i++)
            coeff[i+2] -= (((double) porder) - 1.0) * coeff2[i] / (double) porder;
        
        delete[] coeff1; delete[] coeff2;
    }
    
//     for (i = 0; i <= porder; i++) {
//         k = porder - i;
//         coeff[i] = pow(2.0, (double) porder) * pow(-1.0, (double) k) * nchoosek(porder,k) * nchoosek(0.5*(porder+k-1),porder);
//     }
}

void jacobi(double *coeff, int porder, double alpha, double beta)
{
    //error("jacobi not validated yet.\n");
    
    int i;
    double C1, C2, C3, C4, porder_d = (double) porder;
    
    C1 =  2.0 * porder_d * (porder_d+alpha+beta) * (2.0*porder_d+alpha+beta-2.0);
    C2 =  (2.0*porder_d+alpha+beta-1.0) * (2.0*porder_d+alpha+beta) * (2.0*porder_d+alpha+beta-2.0);
    C3 =  (2.0*porder_d+alpha+beta-1.0) * (alpha*alpha-beta*beta);
    C4 = -2.0 * (porder_d+alpha-1.0) * (porder_d+beta-1) * (2.0*porder_d+alpha+beta);
    
    if (porder == 0)
        coeff[0] = 1.0;
    else if (porder == 1) {
        coeff[0] = 1.0 + 0.5*(alpha+beta);
        coeff[1] = 0.5*(alpha-beta);
    }
    else {
        double *coeff1 = new double[porder];
        double *coeff2 = new double[porder-1];
        
        jacobi(&coeff1[0], porder-1, alpha, beta);
        jacobi(&coeff2[0], porder-2, alpha, beta);
        
        for (i = 0; i <= porder; i++)
            coeff[i] = 0.0;
        for (i = 0; i <= porder-1; i++)
            coeff[i] += C2 * coeff1[i] / C1;
        for (i = 0; i <= porder-1; i++)
            coeff[i+1] += C3 * coeff1[i] / C1;
        for (i = 0; i <= porder-2; i++)
            coeff[i+2] += C4 * coeff2[i] / C1;
        
        delete[] coeff1; delete[] coeff2;
    }
}

void polyder(double *Dcoeff, double *coeff, int porder)
{
    for (int i = 0; i < porder; i++)
        Dcoeff[i] = ((double) (porder - i) ) * coeff[i];
}

void polyval(double *pval, double *coeff, double *x, int porder, int numPoints)
{
    for (int i = 0; i < numPoints; i++) {
        pval[i] = 0.0;
        for (int j = 0; j < porder+1; j++)
            pval[i] += coeff[j] * pow(x[i], (double) (porder-j));
    }
}

void conv(double *output, double *input1, int lenInput1, double *input2, int lenInput2)
{
    int i, j, lenOutput;
    
    if (lenInput1 < 1 || lenInput2 < 1)
        error("Inputs of conv function must be length >= 1.\n");
    
    lenOutput = lenInput1+lenInput2-1;
    for (i = 0; i < lenOutput; i++)
        output[i] = 0.0;
    
    for (i = 0; i < lenInput1; i++)
        for (j = 0; j < lenInput2; j++)
            output[i+j+1+lenOutput-lenInput1-lenInput2] += input1[i]*input2[j];
}

void koornwinder1d(double *f,double *x, int numPoints, int porder, int derivativesFlag)
{
    // f: numPoints / (porder+1) / 2
    // x: numPoints
    
    int i, j;
    
    double *coeff = new double[porder+1];
    double *Dcoeff = new double[porder];
    double *pval = new double[numPoints];
    double *dpval = new double[numPoints];  
    
    // Map from [0,1] to [-1,1]
    for (i = 0; i < numPoints; i++)
        x[i] = 2.0*x[i] - 1.0;
    
    for (i = 0; i < (porder+1); i++) {
        legendre(&coeff[0], i);
        for (j = 0; j <= i; j++)
            coeff[j] = coeff[j]*sqrt(2.0 * ((double) i) + 1.0);
        polyder(&Dcoeff[0],&coeff[0], i);
        polyval(&pval[0], &coeff[0], &x[0], i, numPoints);
        polyval(&dpval[0], &Dcoeff[0], &x[0], i-1, numPoints);
        for (j = 0; j < numPoints; j++) {
            f[i*numPoints + j] = pval[j];
            if (derivativesFlag == 1)
                f[1*(porder+1)*numPoints + i*numPoints + j] = dpval[j];
        }
    }
    
    if (derivativesFlag == 1) {
        for (i = 0; i < (porder+1)*numPoints; i++)
            f[1*(porder+1)*numPoints + i] = 2.0*f[1*(porder+1)*numPoints + i];
    }
    
    // Map from [-1,1] to [0,1]
    for (i = 0; i < numPoints; i++)
        x[i] = 0.5 * (x[i] + 1.0);
    
    delete[] coeff; delete[] Dcoeff;
    delete[] pval; delete[] dpval;
}

void koornwinder2d(double *f,double *x, int numPoints, int porder, int derivativesFlag)
{
    // f: numPoints / numPoly / nd+1
    // x: numPoints / nd
    
    int i, j, k, nd = 2, numPoly = ((porder+1)*(porder+2)) / 2;
    double fc;
    
    double *xc = new double[numPoints*nd];
    double *e = new double[numPoints*nd];
    double *de1 = new double[numPoints*nd];
    double *coeff_p = new double[porder+1];
    double *coeff_q = new double[porder+1];
    double *tmp_coeff_q = new double[porder+1];
    double *Dcoeff_p = new double[porder];
    double *Dcoeff_q = new double[porder];
    double *pval = new double[numPoints];
    double *qval = new double[numPoints];
    double *dpval = new double[numPoints];
    double *dqval = new double[numPoints];
    int *pq = new int[numPoly*nd];
    
    double convVector[] = {-0.5,0.5};
    
    // Map from [0,0]-[1,0]-[0,1] triangle to [-1,-1]-[1,-1]-[-1,1] triangle
    for (i = 0; i < numPoints*nd; i++)
        x[i] = 2.0*x[i] - 1.0;
    
    pascalindex2d(&pq[0], numPoly);
    
    for (j = 0; j < numPoints; j++) {
        xc[0*numPoints+j] = x[0*numPoints+j];
        xc[1*numPoints+j] = min( 0.99999999, x[1*numPoints+j]);     // To avoid singularity
    }
    // xc: numPoints / nd
    
    for (j = 0; j < numPoints; j++) {
        e[0*numPoints+j] = 2.0 * (1.0+xc[0*numPoints+j]) / (1.0-xc[1*numPoints+j]) - 1.0;
        e[1*numPoints+j] = xc[1*numPoints+j];
        if (x[1*numPoints+j] == 1.0) {
            e[0*numPoints+j] = -1.0;
            e[1*numPoints+j] =  1.0;
        }
    }
    // e: numPoints / nd
    
    // Compute f:
    for (i = 0; i < numPoly; i++) {
        jacobi(&coeff_p[0], pq[0*numPoly+i], 0.0, 0.0);
        jacobi(&coeff_q[0], pq[1*numPoly+i], 2.0*((double) pq[0*numPoly+i])+1.0, 0.0);
        
        for (j = 0; j < pq[0*numPoly+i]; j++) {
            conv(&tmp_coeff_q[0], &convVector[0], 2, &coeff_q[0], pq[1*numPoly+i]+1+j);
            for (k = 0; k < pq[1*numPoly+i]+1+j+1; k++)
                coeff_q[k] = tmp_coeff_q[k];
        }
        
        polyval(&pval[0], &coeff_p[0], &e[0*numPoints], pq[0*numPoly+i], numPoints);
        polyval(&qval[0], &coeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i], numPoints);
        
        // Normalization factor to ensure integration to one    
        fc = sqrt(2.0*(2.0*pq[0*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+1.0));
        
        for (j = 0; j < numPoints; j++)
            f[i*numPoints + j] = fc*pval[j]*qval[j];
    }
    
    // Compute derivatives of f (if necessary):
    if (derivativesFlag == 1) {
        // Use displaced coordinate for derivative evaluation
        for (j = 0; j < numPoints; j++) {
            e[0*numPoints+j] = 2.0 * (1.0+xc[0*numPoints+j]) / (1.0-xc[1*numPoints+j]) - 1.0;
            e[1*numPoints+j] = xc[1*numPoints+j];
            de1[0*numPoints+j] = 2.0 / (1.0-xc[1*numPoints+j]);
            de1[1*numPoints+j] = 2.0 * (1.0+xc[0*numPoints+j]) / ((1.0-xc[1*numPoints+j])*(1.0-xc[1*numPoints+j]));
        }
        
        for (i = 0; i < numPoly; i++) {
            jacobi(&coeff_p[0], pq[0*numPoly+i], 0.0, 0.0);
            jacobi(&coeff_q[0], pq[1*numPoly+i], 2.0*((double) pq[0*numPoly+i])+1.0, 0.0);
            for (j = 0; j < pq[0*numPoly+i]; j++) {
                conv(&tmp_coeff_q[0], &convVector[0], 2, &coeff_q[0], pq[1*numPoly+i]+1+j);
                for (k = 0; k < pq[1*numPoly+i]+1+j+1; k++)
                    coeff_q[k] = tmp_coeff_q[k];
            }

            polyder(&Dcoeff_p[0],&coeff_p[0], pq[0*numPoly+i]);
            polyder(&Dcoeff_q[0],&coeff_q[0], pq[0*numPoly+i]+pq[1*numPoly+i]);

            polyval(&pval[0], &coeff_p[0], &e[0*numPoints], pq[0*numPoly+i], numPoints);
            polyval(&dpval[0], &Dcoeff_p[0], &e[0*numPoints], pq[0*numPoly+i]-1, numPoints);
            
            polyval(&qval[0], &coeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i], numPoints);
            polyval(&dqval[0], &Dcoeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i]-1, numPoints);

            // Normalization factor to ensure integration to one    
            fc = sqrt(2.0*(2.0*pq[0*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+1.0));

            for (j = 0; j < numPoints; j++) {
                f[1*numPoly*numPoints + i*numPoints + j] = fc * dpval[j]*qval[j]*de1[0*numPoints+j];
                f[2*numPoly*numPoints + i*numPoints + j] = fc * (dpval[j]*qval[j]*de1[1*numPoints+j] + pval[j]*dqval[j]);
            }
        }
        
        for (j = 0; j < nd*numPoly*numPoints; j++)
            f[numPoly*numPoints+j] = 2.0 * f[numPoly*numPoints+j];
    }
    
    // Map from [-1,-1]-[1,-1]-[-1,1] triangle to [0,0]-[1,0]-[0,1] triangle
    for (i = 0; i < nd*numPoints; i++)
        x[i] = 0.5 * (x[i] + 1.0);
    
    delete[] xc; delete[] e; delete[] de1;
    delete[] coeff_p; delete[] coeff_q; delete[] tmp_coeff_q;
    delete[] Dcoeff_p; delete[] Dcoeff_q;
    delete[] pval; delete[] qval;
    delete[] dpval; delete[] dqval;
    delete[] pq;
}

void koornwinder3d(double *f,double *x, int numPoints, int porder, int derivativesFlag)
{
    // f: numPoints / numPoly / nd+1
    // x: numPoints / nd
    
    int i, j, k, nd = 3, numPoly = ((porder+1)*(porder+2)*(porder+3)) / 6;
    double fc;
    
    double *xc = new double[numPoints*nd];
    double *e = new double[numPoints*nd];
    double *de1 = new double[numPoints*nd];
    double *de2 = new double[numPoints*nd];
    double *coeff_p = new double[porder+1];
    double *coeff_q = new double[porder+1];
    double *coeff_r = new double[porder+1];
    double *tmp_coeff_q = new double[porder+1];
    double *tmp_coeff_r = new double[porder+1];
    double *Dcoeff_p = new double[porder];
    double *Dcoeff_q = new double[porder];
    double *Dcoeff_r = new double[porder];
    double *pval = new double[numPoints];
    double *qval = new double[numPoints];
    double *rval = new double[numPoints];
    double *dpval = new double[numPoints];
    double *dqval = new double[numPoints];
    double *drval = new double[numPoints];
    int *pq = new int[numPoly*nd];
    
    double convVector[2];
    convVector[0] = -0.5; convVector[1] = 0.5;
    
    // Map from [0,0]-[1,0]-[0,1]-[1,1] tet to [-1,-1]-[1,-1]-[-1,1],[1,1] tet
    for (i = 0; i < numPoints*nd; i++)
        x[i] = 2.0*x[i] - 1.0;
    
    pascalindex3d(&pq[0], numPoly);
    
    for (j = 0; j < numPoints; j++) {
        if (x[1*numPoints+j] + x[2*numPoints+j] != 0.0)
            e[0*numPoints+j] = -2.0 * (1.0+x[0*numPoints+j]) / (x[1*numPoints+j]+x[2*numPoints+j]) - 1.0;
        else
            e[0*numPoints+j] = -1.0;
        if (x[2*numPoints+j] != 1.0)
            e[1*numPoints+j] = 2.0 * (1.0+x[1*numPoints+j]) / (1.0-x[2*numPoints+j]) - 1.0;
        else
            e[1*numPoints+j] = -1.0;
        e[2*numPoints+j] = x[2*numPoints+j];
    }
    // e: numPoints / nd

    // Compute f:
    for (i = 0; i < numPoly; i++) {
        jacobi(&coeff_p[0], pq[0*numPoly+i], 0.0, 0.0);
        jacobi(&coeff_q[0], pq[1*numPoly+i], 2.0*((double) pq[0*numPoly+i])+1.0, 0.0);
        jacobi(&coeff_r[0], pq[2*numPoly+i], 2.0*((double) pq[0*numPoly+i])+2.0*((double) pq[1*numPoly+i])+2.0, 0.0);
        for (j = 0; j < pq[0*numPoly+i]; j++) {
            conv(&tmp_coeff_q[0], &convVector[0], 2, &coeff_q[0], pq[1*numPoly+i]+1+j);
            for (k = 0; k < pq[1*numPoly+i]+1+j+1; k++)
                coeff_q[k] = tmp_coeff_q[k];
        }
        for (j = 0; j < pq[0*numPoly+i]+pq[1*numPoly+i]; j++) {
            conv(&tmp_coeff_r[0], &convVector[0], 2, &coeff_r[0], pq[2*numPoly+i]+1+j);
            for (k = 0; k < pq[2*numPoly+i]+1+j+1; k++)
                coeff_r[k] = tmp_coeff_r[k];
        }
        
        polyval(&pval[0], &coeff_p[0], &e[0*numPoints], pq[0*numPoly+i], numPoints);
        polyval(&qval[0], &coeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i], numPoints);
        polyval(&rval[0], &coeff_r[0], &e[2*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i], numPoints);
        
        // Normalization factor to ensure integration to one    
        fc = sqrt(4.0*(2.0*pq[0*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i]+2.0));
        
        for (j = 0; j < numPoints; j++)
            f[i*numPoints + j] = fc*pval[j]*qval[j]*rval[j];
    }
    
    // Compute derivatives of f (if necessary):
    if (derivativesFlag == 1) {
        // Use displaced coordinate for derivative evaluation
        for (j = 0; j < numPoints; j++) {
            xc[0*numPoints+j] = x[0*numPoints+j];
            xc[1*numPoints+j] = x[1*numPoints+j];
            xc[2*numPoints+j] = x[2*numPoints+j];
            if (x[1*numPoints+j] + x[2*numPoints+j] == 0.0)
                xc[2*numPoints+j] = -1.0e-8 - x[1*numPoints+j];
            if (x[2*numPoints+j] == 1.0)
                xc[2*numPoints+j] = 0.99999999;
        }

        for (j = 0; j < numPoints; j++) {
            e[0*numPoints+j] = -2.0 * (1.0+xc[0*numPoints+j]) / (xc[1*numPoints+j]+xc[2*numPoints+j]) - 1.0;
            e[1*numPoints+j] =  2.0 * (1.0+xc[1*numPoints+j]) / (1.0-xc[2*numPoints+j]) - 1.0;
            e[2*numPoints+j] =  xc[2*numPoints+j];
            de1[0*numPoints+j] = - 2.0 / (xc[1*numPoints+j]+xc[2*numPoints+j]);
            de1[1*numPoints+j] = 2.0 * (1.0+xc[0*numPoints+j]) / ((xc[1*numPoints+j]+xc[2*numPoints+j])*(xc[1*numPoints+j]+xc[2*numPoints+j]));
            de1[2*numPoints+j] = 2.0 * (1.0+xc[0*numPoints+j]) / ((xc[1*numPoints+j]+xc[2*numPoints+j])*(xc[1*numPoints+j]+xc[2*numPoints+j]));
            de2[0*numPoints+j] = 0.0 * xc[0*numPoints+j];
            de2[1*numPoints+j] = 2.0 / (1.0-xc[2*numPoints+j]);
            de2[2*numPoints+j] = 2.0 * (1.0+xc[1*numPoints+j]) / ((1.0-xc[2*numPoints+j])*(1.0-xc[2*numPoints+j]));
        }
        
        for (i = 0; i < numPoly; i++) {
            jacobi(&coeff_p[0], pq[0*numPoly+i], 0.0, 0.0);
            jacobi(&coeff_q[0], pq[1*numPoly+i], 2.0*((double) pq[0*numPoly+i])+1.0, 0.0);
            jacobi(&coeff_r[0], pq[2*numPoly+i], 2.0*((double) pq[0*numPoly+i])+2.0*((double) pq[1*numPoly+i])+2.0, 0.0);
            for (j = 0; j < pq[0*numPoly+i]; j++) {
                conv(&tmp_coeff_q[0], &convVector[0], 2, &coeff_q[0], pq[1*numPoly+i]+1+j);
                for (k = 0; k < pq[1*numPoly+i]+1+j+1; k++)
                    coeff_q[k] = tmp_coeff_q[k];
            }
            for (j = 0; j < pq[0*numPoly+i]+pq[1*numPoly+i]; j++) {
                conv(&tmp_coeff_r[0], &convVector[0], 2, &coeff_r[0], pq[2*numPoly+i]+1+j);
                for (k = 0; k < pq[2*numPoly+i]+1+j+1; k++)
                    coeff_r[k] = tmp_coeff_r[k];
            }

            polyder(&Dcoeff_p[0],&coeff_p[0], pq[0*numPoly+i]);
            polyder(&Dcoeff_q[0],&coeff_q[0], pq[0*numPoly+i]+pq[1*numPoly+i]);
            polyder(&Dcoeff_r[0],&coeff_r[0], pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i]);

            polyval(&pval[0], &coeff_p[0], &e[0*numPoints], pq[0*numPoly+i], numPoints);
            polyval(&dpval[0], &Dcoeff_p[0], &e[0*numPoints], pq[0*numPoly+i]-1, numPoints);

            polyval(&qval[0], &coeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i], numPoints);
            polyval(&dqval[0], &Dcoeff_q[0], &e[1*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i]-1, numPoints);
            
            polyval(&rval[0], &coeff_r[0], &e[2*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i], numPoints);
            polyval(&drval[0], &Dcoeff_r[0], &e[2*numPoints], pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i]-1, numPoints);

            // Normalization factor to ensure integration to one    
            fc = sqrt(4.0*(2.0*pq[0*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+1.0)*(pq[0*numPoly+i]+pq[1*numPoly+i]+pq[2*numPoly+i]+2.0));
            
            for (j = 0; j < numPoints; j++) {
                f[1*numPoly*numPoints + i*numPoints + j] = fc * dpval[j]*qval[j]*rval[j]*de1[0*numPoints+j];
                f[2*numPoly*numPoints + i*numPoints + j] = fc * (dpval[j]*qval[j]*rval[j]*de1[1*numPoints+j] + pval[j]*dqval[j]*rval[j]*de2[1*numPoints+j]);
                f[3*numPoly*numPoints + i*numPoints + j] = fc * (dpval[j]*qval[j]*rval[j]*de1[2*numPoints+j] + pval[j]*dqval[j]*rval[j]*de2[2*numPoints+j] + pval[j]*qval[j]*drval[j]);
            }
        }

        for (j = 0; j < nd*numPoly*numPoints; j++)
            f[numPoly*numPoints+j] = 2.0 * f[numPoly*numPoints+j];
    }
    
    // Map from [-1,-1]-[1,-1]-[-1,1],[1,1] tet to [0,0]-[1,0]-[0,1]-[1,1] tet
    for (i = 0; i < nd*numPoints; i++)
        x[i] = 0.5 * (x[i] + 1.0);
    
    delete[] xc; delete[] e; delete[] de1; delete[] de2;
    delete[] coeff_p; delete[] coeff_q; delete[] coeff_r;
    delete[] tmp_coeff_q; delete[] tmp_coeff_r;
    delete[] Dcoeff_p; delete[] Dcoeff_q; delete[] Dcoeff_r;
    delete[] pval; delete[] qval; delete[] rval;
    delete[] dpval; delete[] dqval; delete[] drval;
    delete[] pq;
}

void koornwinder(double *f, double *x, int numPoints, int porder, int nd, int derivativesFlag)
{
    if (nd == 1)
        koornwinder1d(&f[0], &x[0], numPoints, porder, derivativesFlag);
    else if (nd == 2)
        koornwinder2d(&f[0], &x[0], numPoints, porder, derivativesFlag);
    else if (nd == 3)
        koornwinder3d(&f[0], &x[0], numPoints, porder, derivativesFlag);
    else
        error("Number of dimensions not implemented.\n");
}

void tensorproduct(double *f, double *x, int npoints, int porder, int nd, int derivativesFlag)
{
    // f: npoints / npv / nd+1
    // x: npoints / nd
    
    int i, j, k, l;
    
    if (nd == 1) {
        koornwinder1d(&f[0],&x[0],npoints,porder,derivativesFlag);
    }
    else if (nd == 2) {
        double *g1 = new double[2*npoints*(porder+1)];
        double *g2 = new double[2*npoints*(porder+1)];
        double *gx = &g1[npoints*(porder+1)];
        double *gy = &g2[npoints*(porder+1)];
        koornwinder1d(&g1[0],&x[0*npoints],npoints,porder,derivativesFlag);     // Legendre basis in x direction
        koornwinder1d(&g2[0],&x[1*npoints],npoints,porder,derivativesFlag);     // Legendre basis in y direction
        
        // Perform tensor product to obtain the shape functions and their derivatives on the unit square
        for (i = 0; i < npoints; i++) {
            for (j = 0; j < (porder+1); j++) {
                for (k = 0; k < (porder+1); k++) {
                    f[0*(porder+1)*(porder+1)*npoints + (k*(porder+1)+j)*npoints + i] =  g1[j*npoints+i]*g2[k*npoints+i];
                    if (derivativesFlag == 1) {
                        f[1*(porder+1)*(porder+1)*npoints + (k*(porder+1)+j)*npoints + i] = gx[j*npoints+i]*g2[k*npoints+i];
                        f[2*(porder+1)*(porder+1)*npoints + (k*(porder+1)+j)*npoints + i] = g1[j*npoints+i]*gy[k*npoints+i];
                    }
                }
            }
        }
        
        delete[] g1; delete[] g2;
    }
    else if (nd == 3) {
        double *g1 = new double[2*npoints*(porder+1)];
        double *g2 = new double[2*npoints*(porder+1)];
        double *g3 = new double[2*npoints*(porder+1)];
        double *gx = &g1[npoints*(porder+1)];
        double *gy = &g2[npoints*(porder+1)];
        double *gz = &g3[npoints*(porder+1)];
        koornwinder1d(&g1[0],&x[0*npoints],npoints,porder,derivativesFlag);     // Legendre basis in x direction
        koornwinder1d(&g2[0],&x[1*npoints],npoints,porder,derivativesFlag);     // Legendre basis in y direction
        koornwinder1d(&g3[0],&x[2*npoints],npoints,porder,derivativesFlag);     // Legendre basis in z direction
        
        // Perform tensor product to obtain the shape functions and their derivatives on the unit cube
        for (i = 0; i < npoints; i++) {
            for (j = 0; j < (porder+1); j++) {
                for (k = 0; k < (porder+1); k++) {
                    for (l = 0; l < (porder+1); l++) {
                        f[0*(porder+1)*(porder+1)*(porder+1)*npoints + (l*(porder+1)*(porder+1)+k*(porder+1)+j)*npoints + i] =  g1[j*npoints+i]*g2[k*npoints+i]*g3[l*npoints+i];
                        if (derivativesFlag == 1) {
                            f[1*(porder+1)*(porder+1)*(porder+1)*npoints + (l*(porder+1)*(porder+1)+k*(porder+1)+j)*npoints + i] = gx[j*npoints+i]*g2[k*npoints+i]*g3[l*npoints+i];
                            f[2*(porder+1)*(porder+1)*(porder+1)*npoints + (l*(porder+1)*(porder+1)+k*(porder+1)+j)*npoints + i] = g1[j*npoints+i]*gy[k*npoints+i]*g3[l*npoints+i];
                            f[3*(porder+1)*(porder+1)*(porder+1)*npoints + (l*(porder+1)*(porder+1)+k*(porder+1)+j)*npoints + i] = g1[j*npoints+i]*g2[k*npoints+i]*gz[l*npoints+i];
                        }
                    }
                }
            }
        }
        
        delete[] g1; delete[] g2; delete[] g3;
    }
}

void mkshape(vector<double> &shap, vector<double> &plocal, vector<double> &pts, int npoints, int elemtype, int porder, int nd, int numNodes)
{
    // porder: Polynomial order
    // plocal: Node positions. numNodes / nd
    // pts: Points to evaluate shape fucntions and derivatives. npoints / nd
    // shap: shape function and derivatives. numNodes / npoints / nd+1
    
    int i, j, k, info, inc = 1;
    int nd1 = nd + 1;
    int lwork = numNodes;
    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    int *ipiv = new int[numNodes];
    double *work = new double[lwork];
    double *A = new double[numNodes*numNodes];
    double *nf = new double[npoints*numNodes*nd1];
    
    //shap_p.resize(numNodes*(nd+1)*npoints);
    //double *shap = &shap_p[0];
    
    if (elemtype == 0) {
        koornwinder(&nf[0], &pts[0], npoints, porder, nd, 1);           // Orthogonal shape functions
        koornwinder(&A[0], &plocal[0], numNodes, porder, nd, 0);         // Vandermonde matrix
    }
    else if (elemtype == 1) {
        tensorproduct(&nf[0], &pts[0], npoints, porder, nd, 1);         // Orthogonal shape functions
        tensorproduct(&A[0], &plocal[0], numNodes, porder, nd, 0);      // Vandermonde matrix
    }
    else
        error("Element type not implemented.\n");
    
//     if (save==1) {
//       writearray2file("plocal.bin", plocal.data(), numNodes*nd);
//       writearray2file("pts.bin", pts.data(), npoints*nd);
//       writearray2file("nf.bin", nf, npoints*numNodes*nd1);
//       writearray2file("A.bin", A, numNodes*numNodes);
//     }
    
    // Divide orthogonal shape functions by the Vandermonde matrix to obtain nodal shape functions:
    DGETRF(&numNodes, &numNodes, &A[0], &numNodes, &ipiv[0], &info);
    double *Ainv = &A[0];
    DGETRI(&numNodes, &Ainv[0], &numNodes, &ipiv[0], &work[0], &lwork, &info);
    for (i = 0; i < nd1; i++) {
        DGEMM(&chn, &chn, &npoints, &numNodes, &numNodes, &one, &nf[i*numNodes*npoints], &npoints, &Ainv[0],
              &numNodes, &zero, &shap[i*numNodes*npoints], &npoints);
        // nf: npoints / numNodes / nd+1
        // Ainv: numNodes / numNodes
        // shap: npoints / numNodes / nd+1
    }
    
//     if (save==1) {
//       writearray2file("shapt.bin", shap.data(), npoints*numNodes*nd1);
//       writearray2file("Ainv.bin", Ainv, numNodes*numNodes);
//     }
    
    // Permute shap: "npoints / numNodes / nd+1" -> "numNodes / npoints / nd+1"
    double *shap_tmp = new double[numNodes*npoints*(nd+1)];
    for (i = 0; i < nd1; i++)
        for (j = 0; j < numNodes; j++)
            for (k = 0; k < npoints; k++)
                shap_tmp[i*npoints*numNodes+k*numNodes+j] = shap[i*numNodes*npoints+j*npoints+k];
    for (i = 0; i < numNodes*npoints*nd1; i++)
        shap[i] = shap_tmp[i];
    delete[] shap_tmp;
    
//     if (save==1) {
//       writearray2file("shap.bin", shap.data(), npoints*numNodes*nd1);
//     }
    
    delete[] ipiv; delete[] work;
    delete[] A; delete[] nf;
}

void localbasis(double *phielem, double *phiface, const double *plocvl, const double *plocfc,
                int dim, int elemtype, int nne, int nnf)
{
    int i;

    if (dim == 1) {  // 1D line element
        for (i = 0; i < nnf; ++i) {
            phiface[i] = 1.0;  // scalar constant
        }

        for (i = 0; i < nne; ++i) {
            double xi = plocvl[i];  // since plocvl(:,1)
            phielem[i + 0 * nne] = 1.0 - xi;  // phielem(:,1)
            phielem[i + 1 * nne] = xi;        // phielem(:,2)
        }
    }
    else if (dim == 2 && elemtype == 0) {  // triangle
        for (i = 0; i < nnf; ++i) {
            double xi = plocfc[i];  // plocfc(:,1)
            phiface[i + 0 * nnf] = 1.0 - xi;  // phiface(:,1)
            phiface[i + 1 * nnf] = xi;        // phiface(:,2)
        }

        for (i = 0; i < nne; ++i) {
            double xi  = plocvl[i + 0 * nne];  // plocvl(:,1)
            double eta = plocvl[i + 1 * nne];  // plocvl(:,2)
            phielem[i + 0 * nne] = 1.0 - xi - eta;
            phielem[i + 1 * nne] = xi;
            phielem[i + 2 * nne] = eta;
        }
    }
    else if (dim == 2 && elemtype == 1) {  // quadrilateral
        for (i = 0; i < nnf; ++i) {
            double xi = plocfc[i + 0 * nnf];
            phiface[i + 0 * nnf] = 1.0 - xi;
            phiface[i + 1 * nnf] = xi;
        }

        for (i = 0; i < nne; ++i) {
            double xi  = plocvl[i + 0 * nne];
            double eta = plocvl[i + 1 * nne];
            phielem[i + 0 * nne] = (1.0 - xi) * (1.0 - eta);
            phielem[i + 1 * nne] =  xi        * (1.0 - eta);
            phielem[i + 2 * nne] =  xi        * eta;
            phielem[i + 3 * nne] = (1.0 - xi) * eta;
        }
    }
    else if (dim == 3 && elemtype == 0) {  // tetrahedron
        for (i = 0; i < nnf; ++i) {
            double xi  = plocfc[i + 0 * nnf];
            double eta = plocfc[i + 1 * nnf];
            phiface[i + 0 * nnf] = 1.0 - xi - eta;
            phiface[i + 1 * nnf] = xi;
            phiface[i + 2 * nnf] = eta;
        }

        for (i = 0; i < nne; ++i) {
            double xi   = plocvl[i + 0 * nne];
            double eta  = plocvl[i + 1 * nne];
            double zeta = plocvl[i + 2 * nne];
            phielem[i + 0 * nne] = 1.0 - xi - eta - zeta;
            phielem[i + 1 * nne] = xi;
            phielem[i + 2 * nne] = eta;
            phielem[i + 3 * nne] = zeta;
        }
    }
    else if (dim == 3 && elemtype == 1) {  // hexahedron
        for (i = 0; i < nnf; ++i) {
            double xi  = plocfc[i + 0 * nnf];
            double eta = plocfc[i + 1 * nnf];
            phiface[i + 0 * nnf] = (1.0 - xi) * (1.0 - eta);
            phiface[i + 1 * nnf] =  xi        * (1.0 - eta);
            phiface[i + 2 * nnf] =  xi        * eta;
            phiface[i + 3 * nnf] = (1.0 - xi) * eta;
        }

        for (i = 0; i < nne; ++i) {
            double xi   = plocvl[i + 0 * nne];
            double eta  = plocvl[i + 1 * nne];
            double zeta = plocvl[i + 2 * nne];
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

int permindex(vector<int>& permind, const double* plocfc, int npf, int dim, int elemtype) 
{
    int ncols_out = 1;
    if (dim == 1) {         
        permind.resize(1);
        permind[0] = 0;
    } 
    else if (dim == 2) {
        permind.resize(npf);        
        for (int i = 0; i < npf; ++i)
            permind[i] = npf - i - 1;
    } 
    else if (dim == 3 && elemtype == 0) {
        ncols_out = 3;
        permind.resize(npf * 3);
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
        permind.resize(npf * 4);
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

Master initializeMaster(PDE& pde, Mesh& mesh)
{    
    Master master;
    
    std::string fn1 = make_path(pde.exasimpath, "/text2code/text2code/masternodes.bin");
    std::string fn2 = make_path(pde.exasimpath, "/text2code/text2code/gaussnodes.bin");
    
    masternodes(master.xpe, master.telem, master.xpf, master.tface, master.perm, pde.porder, mesh.dim, mesh.elemtype, fn1);     
    gaussnodes(master.gpe, master.gwe, pde.pgauss, mesh.dim, mesh.elemtype, fn2); 
        
    if (mesh.dim>1) gaussnodes(master.gpf, master.gwf, pde.pgauss, mesh.dim-1, mesh.elemtype, fn2);         
    else {
      master.gpf.resize(1); master.gpf[0] = 0.0;
      master.gwf.resize(1); master.gwf[0] = 1.0;
    }
           
    master.nd = mesh.dim;
    master.porder = pde.porder;
    master.pgauss = pde.pgauss;
    master.elemtype = mesh.elemtype;
    master.nodetype = pde.nodetype;
    master.nfe = mesh.nfe;
    master.nve = mesh.nve;
    master.nvf = mesh.nvf;
    master.npe = master.xpe.size()/mesh.dim;
    master.npf = (mesh.dim > 1) ? master.xpf.size()/(mesh.dim-1) : 1;      
    master.nge = master.gwe.size();
    master.ngf = master.gwf.size();
    
//     print2iarray(master.telem.data(), master.telem.size()/mesh.nve, mesh.nve);
//     print2iarray(master.tface.data(), master.tface.size()/mesh.nvf, mesh.nvf);
//     print2iarray(master.perm.data(), master.npf, mesh.nfe);
//     
//     print2darray(master.xpe.data(), master.npe, mesh.dim);
//     print2darray(master.xpf.data(), master.npf, max(mesh.dim-1,1));
//     print2darray(gpe.data(), nge, mesh.dim);
//     print2darray(gwe.data(), nge, 1);
//     print2darray(gpf.data(), ngf, mesh.dim-1);
//     print2darray(gwf.data(), ngf, 1);
    
    master.shapeg.resize(master.npe * master.nge * (master.nd+1));
    master.shapegt.resize(master.npe * master.nge * (master.nd+1));
    master.shapegw.resize(master.npe * master.nge * (master.nd+1));   
    master.shapen.resize(master.npe * master.npe * (master.nd+1));
    master.shapent.resize(master.npe * master.npe * (master.nd+1));
        
    mkshape(master.shapeg, master.xpe, master.gpe, master.nge, master.elemtype, master.porder, master.nd, master.npe);
    mkshape(master.shapen, master.xpe, master.xpe, master.npe, master.elemtype, master.porder, master.nd, master.npe);        
    
    for (int i=0; i<=master.nd; i++)
      for (int j=0; j<master.nge; j++)
        for (int k=0; k<master.npe; k++) {
          master.shapegt[j + master.nge*k + master.npe*master.nge*i] = master.shapeg[k + master.npe*j + master.npe*master.nge*i];
          master.shapegw[k + master.npe*j + master.npe*master.nge*i] = master.shapeg[k + master.npe*j + master.npe*master.nge*i]*master.gwe[j];
        }
    
    for (int i=0; i<=master.nd; i++)
      for (int j=0; j<master.npe; j++)
        for (int k=0; k<master.npe; k++)
          master.shapent[j + master.npe*k + master.npe*master.npe*i] = master.shapen[k + master.npe*j + master.npe*master.npe*i];
    
    master.shapfg.resize(master.npf * master.ngf * (master.nd));
    master.shapfgt.resize(master.npf * master.ngf * (master.nd));
    master.shapfgw.resize(master.npf * master.ngf * (master.nd));   
    master.shapfn.resize(master.npf * master.npf * (master.nd));
    master.shapfnt.resize(master.npf * master.npf * (master.nd));
    
    if (master.nd > 1) {
      mkshape(master.shapfg, master.xpf, master.gpf, master.ngf, master.elemtype, master.porder, master.nd-1, master.npf);
      mkshape(master.shapfn, master.xpf, master.xpf, master.npf, master.elemtype, master.porder, master.nd-1, master.npf);      
    } else {
      master.shapfg[0] = 1.0;
      master.shapfn[0] = 1.0;
    }
    
    for (int i=0; i<master.nd; i++)
      for (int j=0; j<master.ngf; j++)
        for (int k=0; k<master.npf; k++) {
          master.shapfgt[j + master.ngf*k + master.npf*master.ngf*i] = master.shapfg[k + master.npf*j + master.npf*master.ngf*i];
          master.shapfgw[k + master.npf*j + master.npf*master.ngf*i] = master.shapfg[k + master.npf*j + master.npf*master.ngf*i]*master.gwf[j];
        }
    
    for (int i=0; i<master.nd; i++)
      for (int j=0; j<master.npf; j++)
        for (int k=0; k<master.npf; k++)
          master.shapfnt[j + master.npf*k + master.npf*master.npf*i] = master.shapfn[k + master.npf*j + master.npf*master.npf*i];
    
    
    vector<double> b;
    vector<int> a, c, d;
    masternodes(master.xp1d, a, b, c, d, pde.porder, 1, mesh.elemtype, fn1); 
    gaussnodes(master.gp1d, master.gw1d, pde.pgauss, 1, mesh.elemtype, fn2); 
        
    master.np1d = master.xp1d.size(); 
    master.ng1d = master.gw1d.size();    

    master.shap1dg.resize(master.np1d * master.ng1d * 2);
    master.shap1dgt.resize(master.np1d * master.ng1d * 2);
    master.shap1dgw.resize(master.np1d * master.ng1d * 2);   
    master.shap1dn.resize(master.np1d * master.np1d * 2);
    master.shap1dnt.resize(master.np1d * master.np1d * 2);
    
    mkshape(master.shap1dg, master.xp1d, master.gp1d, master.ng1d, master.elemtype, master.porder, 1, master.np1d);
    mkshape(master.shap1dn, master.xp1d, master.xp1d, master.np1d, master.elemtype, master.porder, 1, master.np1d);      
    
    for (int i=0; i<2; i++)
      for (int j=0; j<master.ng1d; j++)
        for (int k=0; k<master.np1d; k++) {
          master.shap1dgt[j + master.ng1d*k + master.np1d*master.ng1d*i] = master.shap1dg[k + master.np1d*j + master.np1d*master.ng1d*i];
          master.shap1dgw[k + master.np1d*j + master.np1d*master.ng1d*i] = master.shap1dg[k + master.np1d*j + master.np1d*master.ng1d*i]*master.gw1d[j];
        }
    
    for (int i=0; i<2; i++)
      for (int j=0; j<master.np1d; j++)
        for (int k=0; k<master.np1d; k++)
          master.shap1dnt[j + master.np1d*k + master.np1d*master.np1d*i] = master.shap1dn[k + master.np1d*j + master.np1d*master.np1d*i];
    
    //for (int i=0; i<=master.nd; i++) print2darray(&master.shapeg[master.npe*master.nge*i], master.npe, master.nge);
    
    master.phielem.resize(master.npe*master.nve);
    master.phiface.resize(master.npf*master.nvf);
    localbasis(master.phielem.data(), master.phiface.data(), master.xpe.data(), master.xpf.data(), master.nd, master.elemtype, master.npe, master.npf);    
    master.npermind = permindex(master.permind, master.xpf.data(), master.npf, master.nd, master.elemtype);              
    
    std::cout << "Finished initializing Master.\n";
           
    return master;
}
         
void writemaster(const Master& master, const std::string& filename) 
{
    //std::cout << "Writing master into file..." << std::endl;

    std::ofstream file(filename, std::ios::binary);
    if (!file) error("Failed to open file for writing: " + filename);
        
    std::vector<double> ndims(20, 0.0);
    ndims[0] = master.nd;
    ndims[1] = master.elemtype;
    ndims[2] = master.nodetype;
    ndims[3] = master.porder;
    ndims[4] = master.pgauss;
    ndims[5] = master.npe;
    ndims[6] = master.npf;
    ndims[7] = master.nge;
    ndims[8] = master.ngf;
    ndims[9] = master.np1d;
    ndims[10] = master.ng1d;

    std::vector<double> nsize(22, 0.0);
    nsize[0] = ndims.size();
    nsize[1] = master.shapegt.size();
    nsize[2] = master.shapegw.size();
    nsize[3] = master.shapfgt.size();
    nsize[4] = master.shapfgw.size();
    nsize[5] = master.shapent.size();
    nsize[6] = master.shapen.size();
    nsize[7] = master.shapfnt.size();
    nsize[8] = master.shapfn.size();
    nsize[9] = master.xpe.size();
    nsize[10] = master.gpe.size();
    nsize[11] = master.gwe.size();
    nsize[12] = master.xpf.size();
    nsize[13] = master.gpf.size();
    nsize[14] = master.gwf.size();
    nsize[15] = master.shap1dgt.size();
    nsize[16] = master.shap1dgw.size();
    nsize[17] = master.shap1dnt.size();
    nsize[18] = master.shap1dn.size();
    nsize[19] = master.xp1d.size();
    nsize[20] = master.gp1d.size();
    nsize[21] = master.gw1d.size();

    auto writeVector = [&](const std::vector<double>& vec) {
        file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));
    };

    double count = static_cast<double>(nsize.size());
    file.write(reinterpret_cast<const char*>(&count), sizeof(double));
    writeVector(nsize);
    writeVector(ndims);
    writeVector(master.shapegt);
    writeVector(master.shapegw);
    writeVector(master.shapfgt);
    writeVector(master.shapfgw);
    writeVector(master.shapent);
    writeVector(master.shapen);
    writeVector(master.shapfnt);
    writeVector(master.shapfn);
    writeVector(master.xpe);
    writeVector(master.gpe);
    writeVector(master.gwe);
    writeVector(master.xpf);
    writeVector(master.gpf);
    writeVector(master.gwf);
    writeVector(master.shap1dgt);
    writeVector(master.shap1dgw);
    writeVector(master.shap1dnt);
    writeVector(master.shap1dn);
    writeVector(master.xp1d);
    writeVector(master.gp1d);
    writeVector(master.gw1d);

    file.close();
    
    std::cout << "Finished writing master to " << filename << std::endl;
}

void buildMesh(Mesh& mesh, const PDE& pde, const Master& master)
{
    mesh.npe = master.npe;
    mesh.npf = master.npf;
    
    if (pde.xdgfile == "") {
        mesh.xdg.resize(master.npe*mesh.dim*mesh.ne);
        compute_dgnodes(mesh.xdg.data(), mesh.p.data(), mesh.t.data(), master.phielem.data(), master.npe, mesh.dim, mesh.ne, mesh.nve);      
        std::cout << "Finished compute_dgnodes.\n";
    }      
    
    mesh.f.resize(mesh.nfe * mesh.ne);
    mesh.t2lf.resize(mesh.nvf * mesh.nfe * mesh.ne);
    mesh.localfaces.resize(mesh.nvf * mesh.nfe);        
    mesh.nf = setboundaryfaces(mesh.f.data(), mesh.t2lf.data(), mesh.localfaces.data(), mesh.p.data(), 
               mesh.t.data(), mesh.boundaryExprs, mesh.dim, mesh.elemtype, mesh.ne, mesh.nbndexpr); 
    
    std::cout << "Finished setboundaryfaces.\n";

    //print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
        
    if (mesh.nprdexpr > 0) {        
        setperiodicfaces(mesh.f.data(), mesh.t.data(), mesh.p.data(), mesh.t2lf.data(), 
            mesh.periodicBoundaries1.data(), mesh.periodicBoundaries2.data(), mesh.periodicExprs1, 
            mesh.periodicExprs2, mesh.dim, mesh.elemtype, mesh.np, mesh.ne, mesh.nprdexpr, mesh.nprdcom); 
        std::cout << "Finished setperiodicfaces.\n";
    }
    
    //print2iarray(mesh.f.data(), mesh.nfe, mesh.ne);
    //print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
    
    if (pde.coupledinterface>0) {
        interface_elements(mesh.inte, mesh.intl, mesh.f, mesh.nfe, mesh.ne, pde.coupledinterface);
        std::cout << "Finished interface_elements.\n";
    }
         
    if (pde.xdgfile == "") {      
      project_dgnodes_onto_curved_boundaries(mesh.xdg.data(), mesh.f.data(), master.perm.data(), mesh.curvedBoundaries.data(),
              mesh.curvedBoundaryExprs, mesh.dim, master.porder, master.npe, master.npf, mesh.nfe, mesh.ne);       
      std::cout << "Finished project_dgnodes_onto_curved_boundaries.\n";
//       print2iarray(mesh.t.data(), mesh.nve, mesh.ne);
//       print2darray(mesh.p.data(), mesh.dim, mesh.np);  
//       print2darray(mesh.xdg.data(), master.npe, mesh.dim);  
//       print2darray(master.phielem.data(), master.npe, mesh.nve);  
    }    
        
    std::cout << "Finished building Mesh.\n";
}

#endif
