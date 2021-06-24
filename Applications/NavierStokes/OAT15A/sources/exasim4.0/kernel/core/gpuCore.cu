#ifndef __GPUCORE
#define __GPUCORE

#include <stdio.h>
#include "gpuArrayOperations.cu"
//#include "gpuGEMM.cu"
#include "gpuArrayPermute.cu"
#include "gpuArrayMultiplication.cu"
#include "gpuElemFaceNode.cu"
#include "gpuElemFaceGeom.cu"
#include "gpuApplyGeom.cu"
#include "gpuApplyFhat.cu"

template <typename T>
__global__ void gpuTemplateApplyGivensRotation(T *H, T *s, T *cs, T *sn,  int i)
{        
    T temp;    
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

template <typename T>
void gpuApplyGivensRotation(T *H, T *s, T *cs, T *sn,  int i)
{    
    gpuTemplateApplyGivensRotation<<<1, 1>>>(H, s, cs, sn, i);
}
template void gpuApplyGivensRotation(double*, double*, double*, double*, int);
template void gpuApplyGivensRotation(float*, float*, float*, float*, int);


template <typename T>
__global__ void gpuTemplateBackSolve(T *y, T *H, T *s, int i, int n)
{
    for (int j=i; j>=0; j--)
        y[j] = s[j];    
    
    for (int j=i; j>=0; j--) {
        y[j] =  y[j]/H[j+n*j]; 
        for (int k=j-1; k>=0; k--)
            y[k] = y[k] - H[k+n*j]*y[j]; 
    }
}

template <typename T>
void gpuBackSolve(T *y, T *H, T *s, int i, int n)
{    
    gpuTemplateBackSolve<<<1, 1>>>(y, H, s, i, n);
}
template void gpuBackSolve(double*, double*, double*, int, int);
template void gpuBackSolve(float*, float*, float*, int, int);


#endif

