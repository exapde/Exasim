#ifndef __GEOMNODE
#define __GEOMNODE

#include <math.h>
#include <omp.h>

template <typename T> void cpuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2)
{        
    int nn = (e2-e1)*np;
    int ncu = nc2-nc1;
    int N = nn*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;
        int j = (idx-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        un[idx] = u[k+(j+nc1)*np+e*np*nc];        
    }        
}

template <typename T> void cpuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2)
{
    int nn = (e2-e1)*np;
    int ncu = nc2-nc1;
    int N = nn*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nn;
        int j = (idx-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        u[k+(j+nc1)*np+e*np*nc] = un[idx];        
    }            
}

template <typename T> void cpuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        #pragma omp parallel for
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
            uh[idx] = 0.5*(udg[m1+j*npe+n1*npe*nc]+udg[m2+j*npe+n2*npe*nc]);
        }                        
    }
    else if (opts==1) {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;
            uh[idx] = udg[m1+j*npe+n1*npe*nc];
        }                        
    }
    else if (opts==2) {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)       
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k2 = facecon[2*m+1];
            int m2 = k2%npe;
            int n2 = (k2-m2)/npe;            
            uh[idx] = udg[m2+j*npe+n2*npe*nc];
        }                        
    }
}


template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts)
{    
    int nf = f2-f1;
    int ndf = npf*nf;
    int N = ndf*ncu;
    
    if (opts==0) {
        #pragma omp parallel for        
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
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
            udg[m2+j*npe+n2*npe*nc] = udg[m2+j*npe+n2*npe*nc] + uh[idx];            
        }                        
    }
    else {
        #pragma omp parallel for
        for (int idx = 0; idx<N; idx++)    
        {
            int i = idx%ndf;
            int j = (idx-i)/ndf;
            int m = npf*f1+i;
            int k1 = facecon[2*m];
            int m1 = k1%npe;
            int n1 = (k1-m1)/npe;            
            udg[m1+j*npe+n1*npe*nc] = udg[m1+j*npe+n1*npe*nc] - uh[idx];
        }                        
    }
}

template <typename T> void cpuElemGeom(T *Xx, T *jac, T *Jg, int ne, int ng, int nd)
{        
    T *Jg11, *Jg12, *Jg13, *Jg21, *Jg22, *Jg23, *Jg31, *Jg32, *Jg33;
    T *Xx11, *Xx12, *Xx13, *Xx21, *Xx22, *Xx23, *Xx31, *Xx32, *Xx33;
    int ngv = ng*ne;
     
    if (nd==1) {
        #pragma omp parallel for
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
        #pragma omp parallel for
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
        #pragma omp parallel for
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

template <typename T> void cpuFaceGeom(T *nlg, T *jacg, T *Jg, int nf, int ng, int nd)
{        
    T *Jg11, *Jg12, *Jg21, *Jg22, *Jg31, *Jg32;
    int na = ng*nf;        
    
    if (nd==1) {
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            jacg[i] = 1.0;
            nlg[i] = -1.0;
        }                
    }
    else if (nd==2) {        
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            int j = i+na;
            jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
            nlg[i] = Jg[j]/jacg[i];
            nlg[j] = -Jg[i]/jacg[i];
        }
    }
    else if (nd==3) {
        Jg11 = &Jg[0];
        Jg21 = &Jg[na];
        Jg31 = &Jg[2*na];
        Jg12 = &Jg[3*na];
        Jg22 = &Jg[4*na];
        Jg32 = &Jg[5*na];
        #pragma omp parallel for
        for (int i=0; i<na; i++) {
            int j = i+na;
            int k = i+2*na;
            nlg[i] = Jg21[i]*Jg32[i] - Jg31[i]*Jg22[i];
            nlg[j] = Jg31[i]*Jg12[i] - Jg11[i]*Jg32[i];
            nlg[k] = Jg11[i]*Jg22[i] - Jg21[i]*Jg12[i];
            jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
            nlg[i] = nlg[i]/jacg[i];
            nlg[j] = nlg[j]/jacg[i];
            nlg[k] = nlg[k]/jacg[i];
        }
    }    
}

template <typename T> void cpuApplyJac1(T *sg, T *jac, int nga, int ncu)
{
    int N = nga*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nga;        
        sg[idx] = sg[idx]*jac[i];
    }
}

template <typename T> void cpuApplyJac2(T *sg, T *jac, T *ug, T *su, T *fc_u, int nga, int ncu)
{
    int N = nga*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nga;        
        int n = (idx-i)/nga;
        sg[idx] = (sg[idx] + su[idx] - ug[idx]*fc_u[n])*jac[i];
    }    
}

template <typename T> void cpuApplyXx1(T *sg, T *ug, T *Xx, int nga, int nd, int ncu)
{
    int N = nga*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nga;        
        sg[idx] = ug[idx]*Xx[i];
    }    
}

template <typename T> void cpuApplyXx2(T *sg, T *fg, T *Xx, int nga, int nd, int ncu)
{
    int N = nga*ncu;

    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nga;        
        sg[idx] = 0.0;
        for (int j=0; j<nd; j++)
            sg[idx] = sg[idx] + fg[idx+N*j]*Xx[i+nga*j];
    }        
}

template <typename T> void cpuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd)
{
    int N = nga*ncu;
    #pragma omp parallel for
    for (int idx = 0; idx<N; idx++)
    {
        int i = idx%nga;                
        for (int j=0; j<nd; j++)
            fqg[idx+N*j] = uhg[idx]*nlg[i+nga*j]*jac[i];
    }            
}

template void cpuGetElemNodes(double*, double*, int, int, int, int, int, int);
template void cpuPutElemNodes(double*, double*, int, int, int, int, int, int);
template void cpuGetFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(double*, double*, int*, int, int, int, int, int, int, int);
template void cpuElemGeom(double*, double*, double*, int, int, int);
template void cpuFaceGeom(double*, double*, double*, int, int, int);
template void cpuApplyJac1(double*, double*, int, int);
template void cpuApplyJac2(double*, double*, double*, double*, double*, int, int);
template void cpuApplyXx1(double*, double*, double*, int, int, int);
template void cpuApplyXx2(double*, double*, double*, int, int, int);
template void cpuApplyJacNormal(double*, double*, double*, double*, int, int, int);    

template void cpuGetElemNodes(float*, float*, int, int, int, int, int, int);
template void cpuPutElemNodes(float*, float*, int, int, int, int, int, int);
template void cpuGetFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void cpuPutFaceNodes(float*, float*, int*, int, int, int, int, int, int, int);
template void cpuElemGeom(float*, float*, float*, int, int, int);
template void cpuFaceGeom(float*, float*, float*, int, int, int);
template void cpuApplyJac1(float*, float*, int, int);
template void cpuApplyJac2(float*, float*, float*, float*, float*, int, int);
template void cpuApplyXx1(float*, float*, float*, int, int, int);
template void cpuApplyXx2(float*, float*, float*, int, int, int);
template void cpuApplyJacNormal(float*, float*, float*, float*, int, int, int);    


#endif


