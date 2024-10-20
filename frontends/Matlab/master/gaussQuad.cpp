#ifndef __GAUSSQUAD
#define __GAUSSQUAD

#include "gaussQuadTri.cpp"
#include "gaussQuadTet.cpp"

// Written by: C. Nguyen & P. Fernandez

void gaussQuad1d(double *x, double *w, Int pgauss)
{
    Int i, ng, lwork, info;
    char chv = 'V', chu = 'U';
    
    if (pgauss % 2)     // pgauss is odd
        ng = (pgauss-1)/2 + 1;
    else                // pgauss is even
        ng = pgauss/2 + 1;
    
    if (ng > 1) {
        lwork = 3*ng-1;
        
        double *beta = new double[ng-1];
        double *T = new double[ng*ng];
        double *L = new double[ng];
        double *work = new double[lwork];

        // 3-term recurrence coefficients for Jacobi matrix:
        for (i = 1; i < ng; i++)
            beta[i-1] = 0.5 / sqrt(1.0 - 1.0/((2.0 * (double) i)*(2.0 * (double) i)));

        // Construct upper triangular part of Jacobi matrix T
        for (i = 0; i < ng*ng; i++)
            T[i] = 0.0;
        for (i = 1; i < ng; i++)
            T[(ng+1)*i-1] = beta[i-1];

        DSYEV(&chv, &chu, &ng, &T[0], &ng, &L[0], &work[0], &lwork, &info);

        // Compute nodes in [0,1] (= Legendre points):
        for (i = 0; i < ng; i++)
            x[i] = 0.5*(L[i] + 1.0);

        // Compute weights in [0 1]:
        for (i = 0; i < ng; i++)
            w[i] = T[i*ng]*T[i*ng];

        delete[] beta; delete[] T;
        delete[] L;
        delete[] work;
    }
    else {
        x[0] = 0.5;
        w[0] = 1.0;
    }
}

void gaussQuadTensor(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss, Int dim)
{
    // x: ng / nd
    // w: ng
    
    Int i, j, k, ng1d, ng2d, ng3d;
    
    if (pgauss % 2)     // pgauss is odd
        ng1d = (pgauss-1)/2 + 1;
    else                // pgauss is even
        ng1d = pgauss/2 + 1;
    ng2d = ng1d*ng1d;
    ng3d = ng1d*ng1d*ng1d;
    
    switch (dim) {
        case 1:
            *nq = ng1d;
            break;
        case 2:
            *nq = ng2d;
            break;
        case 3:
            *nq = ng3d;
            break;
        default:
            error("Dimension not implemented.\n");
    }
    
    x_p[0].resize((*nq)*dim);
    w_p[0].resize((*nq));
    
    double *x = &x_p[0][0];
    double *w = &w_p[0][0];
    
    double *x1d = new double[ng1d];
    double *w1d = new double[ng1d];
    
    gaussQuad1d(&x1d[0], &w1d[0], pgauss);
    
    switch (dim) {
        case 1:
            for (i = 0; i < ng1d; i++) {
                x[i] = x1d[i];
                w[i] = w1d[i];
            }
            break;
        case 2:
            for (i = 0; i < ng1d; i++) {
                for (j = 0; j < ng1d; j++) {
                    x[0*ng2d+i*ng1d+j] = x1d[j];
                    x[1*ng2d+i*ng1d+j] = x1d[i];
                    w[i*ng1d+j] = w1d[i]*w1d[j];
                }
            }        
            break;
        case 3:
            for (i = 0; i < ng1d; i++) {
                for (j = 0; j < ng1d; j++) {
                    for (k = 0; k < ng1d; k++) {
                        x[0*ng3d+i*ng2d+j*ng1d+k] = x1d[k];
                        x[1*ng3d+i*ng2d+j*ng1d+k] = x1d[j];
                        x[2*ng3d+i*ng2d+j*ng1d+k] = x1d[i];
                        w[i*ng2d+j*ng1d+k] = w1d[i]*w1d[j]*w1d[k];
                    }
                }
            }   
            break;
        default:
            error("Dimension not implemented.\n");
    }
    
    delete[] x1d; delete[] w1d;
}

void gaussQuad(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss, Int dim, Int elemtype)
{
    switch (dim) {
        case 0:
            *nq = 1;
            x_p[0].resize(1);
            x_p[0][0] = 0.0;
            w_p[0].resize(1);
            w_p[0][0] = 1.0;
            break;
        case 1:
            gaussQuadTensor(x_p, w_p, nq, pgauss, dim);
            break;
        case 2:
            if (elemtype == 0)
                gaussQuadTri(x_p, w_p, nq, pgauss);          // tri
            else if (elemtype == 1)
                gaussQuadTensor(x_p, w_p, nq, pgauss, dim);     // quad
            break;
        case 3:
            if (elemtype == 0)
                gaussQuadTet(x_p, w_p, nq, pgauss);          // tet
            else if (elemtype == 1)
                gaussQuadTensor(x_p, w_p, nq, pgauss, dim);     // hex
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

#endif
