#ifndef __KOORNWINDER
#define __KOORNWINDER

// Written by: C. Nguyen & P. Fernandez

void pascalindex2d(Int *pq, Int numPoly)
{
    Int i, j, counter = 0, done = 0;
    
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

void pascalindex3d(Int *pq, Int numPoly)
{
    Int i, j, k, counter = 0, done = 0;
    
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

void legendre(double *coeff, Int porder)
{
    Int i;
    
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

void jacobi(double *coeff, Int porder, double alpha, double beta)
{
    //error("jacobi not validated yet.\n");
    
    Int i;
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

void polyder(double *Dcoeff, double *coeff, Int porder)
{
    for (Int i = 0; i < porder; i++)
        Dcoeff[i] = ((double) (porder - i) ) * coeff[i];
}

void polyval(double *pval, double *coeff, double *x, Int porder, Int numPoints)
{
    for (Int i = 0; i < numPoints; i++) {
        pval[i] = 0.0;
        for (Int j = 0; j < porder+1; j++)
            pval[i] += coeff[j] * pow(x[i], (double) (porder-j));
    }
}

void conv(double *output, double *input1, Int lenInput1, double *input2, Int lenInput2)
{
    Int i, j, lenOutput;
    
    if (lenInput1 < 1 || lenInput2 < 1)
        error("Inputs of conv function must be length >= 1.\n");
    
    lenOutput = lenInput1+lenInput2-1;
    for (i = 0; i < lenOutput; i++)
        output[i] = 0.0;
    
    for (i = 0; i < lenInput1; i++)
        for (j = 0; j < lenInput2; j++)
            output[i+j+1+lenOutput-lenInput1-lenInput2] += input1[i]*input2[j];
}

void koornwinder1d(double *f,double *x, Int numPoints, Int porder, Int derivativesFlag)
{
    // f: numPoints / (porder+1) / 2
    // x: numPoints
    
    Int i, j;
    
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

void koornwinder2d(double *f,double *x, Int numPoints, Int porder, Int derivativesFlag)
{
    // f: numPoints / numPoly / nd+1
    // x: numPoints / nd
    
    Int i, j, k, nd = 2, numPoly = ((porder+1)*(porder+2)) / 2;
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
    Int *pq = new Int[numPoly*nd];
    
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

void koornwinder3d(double *f,double *x, Int numPoints, Int porder, Int derivativesFlag)
{
    // f: numPoints / numPoly / nd+1
    // x: numPoints / nd
    
    Int i, j, k, nd = 3, numPoly = ((porder+1)*(porder+2)*(porder+3)) / 6;
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
    Int *pq = new Int[numPoly*nd];
    
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

void koornwinder(double *f, double *x, Int numPoints, Int porder, Int nd, Int derivativesFlag)
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

#endif
