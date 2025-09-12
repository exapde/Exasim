#ifndef __GAUSSLOBATTOQUAD
#define __GAUSSLOBATTOQUAD

// Written by: C. Nguyen & P. Fernandez

void gaussLobattoQuad1d(double *x, double *w, Int pgauss)
{
    // x: nq
    // w: nq
    // Ref: http://mathworld.wolfram.com/LobattoQuadrature.html
    
    Int nq;
    if (pgauss % 2)     // pgauss is odd
        nq = (pgauss+3)/2;
    else                // pgauss is even
        nq = (pgauss+4)/2;
    
    switch (nq) {
        case 1:
            error("Gauss Lobatto quadrature not applicable for one quadrature point only...\n");
            break;
        case 2:
            x[0] = 0.0;
            x[1] = 1.0;
            w[0] = 0.5;
            w[1] = 0.5;
            break;
        case 3:
            x[0] = 0.0;
            x[1] = 0.5;
            x[2] = 1.0;
            w[0] = 0.5/3.0;
            w[1] = 2.0/3.0;
            w[2] = 0.5/3.0;
            break;
        case 4:
            x[0] =  0.0;
            x[1] = -0.5 * sqrt(5.0) / 5.0 + 0.5;
            x[2] =  0.5 * sqrt(5.0) / 5.0 + 0.5;
            x[3] =  1.0;
            w[0] = 1.0 / 12.0;
            w[1] = 5.0 / 12.0;
            w[2] = 5.0 / 12.0;
            w[3] = 1.0 / 12.0;
            break;
        case 5:
            x[0] =  0.0;
            x[1] = -0.5*sqrt(21.0) / 7.0 + 0.5;
            x[2] =  0.5;
            x[3] =  0.5*sqrt(21.0) / 7.0 + 0.5;
            x[4] =  1.0;
            w[0] = 1.0 / 20.0;
            w[1] = 49.0 / 180.0;
            w[2] = 16.0 / 45.0;
            w[3] = 49.0 / 180.0;
            w[4] = 1.0 / 20.0;
            break;
        case 6:
            x[0] =  0.0;
            x[1] = -0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[2] = -0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[3] =  0.5 * sqrt((7.0-2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[4] =  0.5 * sqrt((7.0+2.0*sqrt(7.0)) / 21.0) + 0.5;
            x[5] =  1.0;
            w[0] = 1.0 / 30.0;
            w[1] = (14.0-sqrt(7.0)) / 60.0;
            w[2] = (14.0+sqrt(7.0)) / 60.0;
            w[3] = (14.0+sqrt(7.0)) / 60.0;
            w[4] = (14.0-sqrt(7.0)) / 60.0;
            w[5] = 1.0 / 30.0;
            break;
        default:
            error("Gauss-Lobatto quadrature not implemented for pgauss > 9.\n");
    }
}

void gaussLobattoQuadTensor(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss, Int dim)
{
    // x: nq / nd
    // w: nq
    
    Int i, j, k, nq1d, nq2d, nq3d;
    
    if (pgauss % 2)     // pgauss is odd
        nq1d = (pgauss+3)/2;
    else                // pgauss is even
        nq1d = (pgauss+4)/2;
    nq2d = nq1d*nq1d;
    nq3d = nq1d*nq1d*nq1d;
    
    switch (dim) {
        case 1:
            *nq = nq1d;
            break;
        case 2:
            *nq = nq2d;
            break;
        case 3:
            *nq = nq3d;
            break;
        default:
            error("Dimension not implemented.\n");
    }
    
    x_p[0].resize((*nq)*dim);
    w_p[0].resize((*nq));
    
    double *x = &x_p[0][0];
    double *w = &w_p[0][0];
    
    double *x1d = new double[nq1d];
    double *w1d = new double[nq1d];
    
    gaussLobattoQuad1d(&x1d[0], &w1d[0], pgauss);
    
    switch (dim) {
        case 1:
            for (i = 0; i < nq1d; i++) {
                x[i] = x1d[i];
                w[i] = w1d[i];
            }
            break;
        case 2:
            for (i = 0; i < nq1d; i++) {
                for (j = 0; j < nq1d; j++) {
                    x[0*nq2d+i*nq1d+j] = x1d[j];
                    x[1*nq2d+i*nq1d+j] = x1d[i];
                    w[i*nq1d+j] = w1d[i]*w1d[j];
                }
            }        
            break;
        case 3:
            for (i = 0; i < nq1d; i++) {
                for (j = 0; j < nq1d; j++) {
                    for (k = 0; k < nq1d; k++) {
                        x[0*nq3d+i*nq2d+j*nq1d+k] = x1d[k];
                        x[1*nq3d+i*nq2d+j*nq1d+k] = x1d[j];
                        x[2*nq3d+i*nq2d+j*nq1d+k] = x1d[i];
                        w[i*nq2d+j*nq1d+k] = w1d[i]*w1d[j]*w1d[k];
                    }
                }
            }   
            break;
        default:
            error("Dimension not implemented.\n");
    }
    
    delete[] x1d; delete[] w1d;
}

void gaussLobattoQuad(vector<double> *x_p, vector<double> *w_p, Int *nq, Int pgauss, Int dim, Int elemtype)
{
    printf("Gauss-Lobatto quadrature not validated yet.\n");
    
    switch (dim) {
        case 0:
            *nq = 1;
            x_p[0].resize(1);
            x_p[0][0] = 0.0;
            w_p[0].resize(1);
            w_p[0][0] = 1.0;
            break;
        case 1:
            gaussLobattoQuadTensor(x_p, w_p, nq, pgauss, dim);
            break;
        case 2:
            if (elemtype == 0)
                error("Gauss-Lobatto quadrature not implemented for tris.\n");
            else if (elemtype == 1)
                gaussLobattoQuadTensor(x_p, w_p, nq, pgauss, dim);     // quad
            break;
        case 3:
            if (elemtype == 0)
                error("Gauss-Lobatto quadrature not implemented for tets.\n");
            else if (elemtype == 1)
                gaussLobattoQuadTensor(x_p, w_p, nq, pgauss, dim);     // hex
            break;
        default:
            error("Dimension not implemented.\n");
    }
}

#endif
