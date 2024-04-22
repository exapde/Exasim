#ifndef __MKSHAPE
#define __MKSHAPE

#include "koornwinder.cpp"
#include "tensorproduct.cpp"

// Written by: C. Nguyen & P. Fernandez

void mkshape(vector<double> *shap_p, double *plocal, double *pts, Int npoints, Int elemtype, Int porder, Int nd, Int numNodes)
{
    // porder: Polynomial order
    // plocal: Node positions. numNodes / nd
    // pts: Points to evaluate shape fucntions and derivatives. npoints / nd
    // shap: shape function and derivatives. numNodes / npoints / nd+1
    
    Int i, j, k, info, inc = 1;
    Int nd1 = nd + 1;
    Int lwork = numNodes;
    char chn = 'N';
    double one = 1.0, zero = 0.0;
    
    Int *ipiv = new Int[numNodes];
    double *work = new double[lwork];
    double *A = new double[numNodes*numNodes];
    double *nf = new double[npoints*numNodes*nd1];
    
    shap_p[0].resize(numNodes*(nd+1)*npoints);
    double *shap = &shap_p[0][0];
    
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
    
    // Permute shap: "npoints / numNodes / nd+1" -> "numNodes / npoints / nd+1"
    double *shap_tmp = new double[numNodes*npoints*(nd+1)];
    for (i = 0; i < nd1; i++)
        for (j = 0; j < numNodes; j++)
            for (k = 0; k < npoints; k++)
                shap_tmp[i*npoints*numNodes+k*numNodes+j] = shap[i*numNodes*npoints+j*npoints+k];
    for (i = 0; i < numNodes*npoints*nd1; i++)
        shap[i] = shap_tmp[i];
    delete[] shap_tmp;
    
    delete[] ipiv; delete[] work;
    delete[] A; delete[] nf;
}

#endif
