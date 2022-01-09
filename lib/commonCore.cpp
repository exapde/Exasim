#ifndef __COMMONCORE
#define __COMMONCORE

#include <math.h>
#include <algorithm>

void faceindex2(int *in1, int *in2, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2)
{    
    int nf = f2-f1;
    int ndf = npf*ncu;
    int N = ndf*nf;    
    for (int idx = 0; idx<N; idx++)
    {        
        int i = idx%ndf;     // [1, npf*ncu]
        int e = (idx-i)/ndf; // [1, nf]
        int g = i%npf;       // [1, npf]
        int j = (i-g)/npf;   // [1, ncu]                
        int m = g + npf*(f1+e);
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

void faceperm2(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
{
    int N = 0;
    for (int j=0; j<nbf; j++) {
        int f1 = fblks[3*j]-1;
        int f2 = fblks[3*j+1];      
        int nf = f2-f1;             
        faceindex2(&ind1[N], &ind2[N], facecon, npf, ncu, npe, nc, f1, f2);
        indpts[j] = N;
        N = N + npf*ncu*nf;    
    }        
    indpts[nbf] = N;
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

// void mysort(int *a, int n)
// {
//     int b[n];
//     for (int i=0; i<n; i++)
//         b[i] = a[i];
//     std::sort(b,b+n);
//     for (int i=0; i<n; i++)
//         a[i] = b[i];
// }
// 
// void faceperm(int *ind1, int *ind2, int *indpts, int *facecon, int *fblks, int npf, int ncu, int npe, int nc, int nbf)
// {
//     int N = 0;
//     for (int j=0; j<nbf; j++) {
//         int f1 = fblks[3*j]-1;
//         int f2 = fblks[3*j+1];      
//         int nf = f2-f1;
//         int ndf = npf*nf;      
//         int P = N+ndf*ncu;
//         faceindex(&ind1[N], &ind2[N], facecon, npf, ncu, npe, nc, f1, f2);
//         mysort(&ind1[N],ndf*ncu);
//         mysort(&ind2[N],ndf*ncu);
//         indpts[j] = N;
//         N = N + ndf*ncu;    
//     }        
//     indpts[nbf] = N;
// }

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

template <typename T> T cpuArrayGetElementAtIndex(T *y, int n)
{        
    return y[n];        
}
template double cpuArrayGetElementAtIndex(double*, int);
template float cpuArrayGetElementAtIndex(float*, int);

template <typename T> void cpuArraySetValueAtIndex(T *y, T a, int n)
{        
    y[n] = a;        
}
template void cpuArraySetValueAtIndex(double*, double, int);
template void cpuArraySetValueAtIndex(float*, float, int);

template <typename T> void cpuApplyGivensRotation(T *H, T *s, T *cs, T *sn,  int i)
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
template void cpuApplyGivensRotation(double*, double*, double*, double*, int);
template void cpuApplyGivensRotation(float*, float*, float*, float*, int);

template <typename T> void cpuBackSolve(T *y, T *H, T *s, int i, int n)
{
    for (int j=i; j>=0; j--)
        y[j] = s[j];    
    
    for (int j=i; j>=0; j--) {
        y[j] =  y[j]/H[j+n*j]; 
        for (int k=j-1; k>=0; k--)
            y[k] = y[k] - H[k+n*j]*y[j]; 
    }
}
template void cpuBackSolve(double*, double*, double*, int, int);
template void cpuBackSolve(float*, float*, float*, int, int);

#endif






