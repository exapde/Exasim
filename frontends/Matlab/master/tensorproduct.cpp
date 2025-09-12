#ifndef __TENSORPRODUCT
#define __TENSORPRODUCT

// Written by: C. Nguyen & P. Fernandez

void tensorproduct(double *f, double *x, Int npoints, Int porder, Int nd, Int derivativesFlag)
{
    // f: npoints / npv / nd+1
    // x: npoints / nd
    
    Int i, j, k, l;
    
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

#endif
