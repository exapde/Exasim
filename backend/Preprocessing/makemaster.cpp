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

#endif
