#ifndef __GPUGEOM
#define __GPUGEOM

template <typename T>
__global__ void gpuTemplateGetElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%nn; // ne*ncu
        int j = (tid-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        int ind = k+(j+nc1)*np+e*np*nc;
        un[tid] = u[ind];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateGetElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc, int gridDim, int blockDim)
{
    gpuTemplateGetElemNodes<<<gridDim, blockDim>>>(un, u, N, nn, np, e1, nc1, nc);
}
template void gpuTemplateGetElemNodes(double*, double*, int, int, int, int, int, int, int, int);
template void gpuTemplateGetElemNodes(float*, float*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void gpuTemplatePutElemNodes(T *u, T *un, int N, int nn, int np, int e1, int nc1, int nc)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%nn;
        int j = (tid-i)/nn;
        int k = i%np;
        int e = (i-k)/np+e1;        
        int ind = k+(j+nc1)*np+e*np*nc;
        u[ind] = un[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplatePutElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc, int gridDim, int blockDim)
{
    gpuTemplatePutElemNodes<<<gridDim, blockDim>>>(un, u, N, nn, np, e1, nc1, nc);
}
template void gpuTemplatePutElemNodes(double*, double*, int, int, int, int, int, int, int, int);
template void gpuTemplatePutElemNodes(float*, float*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void gpuTemplateGetFaceNodes0(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%ndf;
        int j = (tid-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int k2 = facecon[2*m+1];
        int m1 = k1%npe;
        int m2 = k2%npe;
        int n1 = (k1-m1)/npe;
        int n2 = (k2-m2)/npe;          
        int ind1 = m1+j*npe+n1*npe*nc;
        int ind2 = m2+j*npe+n2*npe*nc;
        uh[tid] = 0.5*(udg[ind1]+udg[ind2]);
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateGetFaceNodes0(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1, int gridDim, int blockDim)
{
    gpuTemplateGetFaceNodes0<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuTemplateGetFaceNodes0(double*, double*, int*, int, int, int, int, int, int, int, int);
template void gpuTemplateGetFaceNodes0(float*, float*, int*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplateGetFaceNodes1(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%ndf;
        int j = (tid-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int m1 = k1%npe;
        int n1 = (k1-m1)/npe;
        int ind1 = m1+j*npe+n1*npe*nc;
        uh[tid] = udg[ind1];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateGetFaceNodes1(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1, int gridDim, int blockDim)
{
    gpuTemplateGetFaceNodes1<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuTemplateGetFaceNodes1(double*, double*, int*, int, int, int, int, int, int, int, int);
template void gpuTemplateGetFaceNodes1(float*, float*, int*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplateGetFaceNodes2(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%ndf;
        int j = (tid-i)/ndf;
        int m = npf*f1+i;
        int k2 = facecon[2*m+1];
        int m2 = k2%npe;
        int n2 = (k2-m2)/npe;             
        int ind2 = m2+j*npe+n2*npe*nc;       
        uh[tid] = udg[ind2];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateGetFaceNodes2(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1, int gridDim, int blockDim)
{
    gpuTemplateGetFaceNodes1<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuTemplateGetFaceNodes2(double*, double*, int*, int, int, int, int, int, int, int, int);
template void gpuTemplateGetFaceNodes2(float*, float*, int*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplatePutFaceNodes0(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
        int i = tid%ndf;
        int j = (tid-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int k2 = facecon[2*m+1];
        int m1 = k1%npe;
        int m2 = k2%npe;
        int n1 = (k1-m1)/npe;
        int n2 = (k2-m2)/npe;          
        int ind1 = m1+j*npe+n1*npe*nc;
        int ind2 = m2+j*npe+n2*npe*nc;
        udg[ind1] = udg[ind1] - uh[tid];
        udg[ind2] = udg[ind2] + uh[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplatePutFaceNodes0(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1, int gridDim, int blockDim)
{
    gpuTemplatePutFaceNodes0<<<gridDim, blockDim>>>(udg, uh, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuTemplatePutFaceNodes0(double*, double*, int*, int, int, int, int, int, int, int, int);
template void gpuTemplatePutFaceNodes0(float*, float*, int*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplatePutFaceNodes1(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N) {
        int i = tid%ndf;
        int j = (tid-i)/ndf;
        int m = npf*f1+i;
        int k1 = facecon[2*m];
        int m1 = k1%npe;
        int n1 = (k1-m1)/npe;
        int ind1 = m1+j*npe+n1*npe*nc;        
        udg[ind1] = udg[ind1] - uh[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplatePutFaceNodes1(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1, int gridDim, int blockDim)
{
    gpuTemplatePutFaceNodes1<<<gridDim, blockDim>>>(udg, uh, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuTemplatePutFaceNodes1(double*, double*, int*, int, int, int, int, int, int, int, int);
template void gpuTemplatePutFaceNodes1(float*, float*, int*, int, int, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplateElemGeom1D(T *jac, T *Xx, T *Jg, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg[tid];
        Xx[tid] = 1.0;
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateElemGeom1D(T *jac, T *Xx, T *Jg, int nga, int gridDim, int blockDim)
{
    gpuTemplateElemGeom1D<<<gridDim, blockDim>>>(jac, Xx, Jg, nga);
}
template void gpuTemplateElemGeom1D(double*, double*, double*, int, int, int);
template void gpuTemplateElemGeom1D(float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22,
                 T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg11[tid]*Jg22[tid] - Jg12[tid]*Jg21[tid];
        Xx11[tid] = Jg22[tid];
        Xx21[tid] = -Jg21[tid];
        Xx12[tid] = -Jg12[tid];
        Xx22[tid] = Jg11[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22,
                 T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga, int gridDim, int blockDim)
{
    gpuTemplateElemGeom2D<<<gridDim, blockDim>>>(jac, Xx11, Xx12, Xx21, Xx22, Jg11, Jg12, Jg21, Jg22, nga);
}
template void gpuTemplateElemGeom2D(double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
template void gpuTemplateElemGeom2D(float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, 
                T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33,
                T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, 
                T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg11[tid]*Jg22[tid]*Jg33[tid] - Jg11[tid]*Jg32[tid]*Jg23[tid] +
                 Jg21[tid]*Jg32[tid]*Jg13[tid] - Jg21[tid]*Jg12[tid]*Jg33[tid] +
                 Jg31[tid]*Jg12[tid]*Jg23[tid] - Jg31[tid]*Jg22[tid]*Jg13[tid];
        Xx11[tid] = Jg22[tid]*Jg33[tid] - Jg23[tid]*Jg32[tid];
        Xx21[tid] = Jg23[tid]*Jg31[tid] - Jg21[tid]*Jg33[tid];
        Xx31[tid] = Jg21[tid]*Jg32[tid] - Jg22[tid]*Jg31[tid];
        Xx12[tid] = Jg13[tid]*Jg32[tid] - Jg12[tid]*Jg33[tid];
        Xx22[tid] = Jg11[tid]*Jg33[tid] - Jg13[tid]*Jg31[tid];
        Xx32[tid] = Jg12[tid]*Jg31[tid] - Jg11[tid]*Jg32[tid];
        Xx13[tid] = Jg12[tid]*Jg23[tid] - Jg13[tid]*Jg22[tid];
        Xx23[tid] = Jg13[tid]*Jg21[tid] - Jg11[tid]*Jg23[tid];
        Xx33[tid] = Jg11[tid]*Jg22[tid] - Jg12[tid]*Jg21[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, 
                T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33,
                T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, 
                T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga, int gridDim, int blockDim)
{
    gpuTemplateElemGeom3D<<<gridDim, blockDim>>>(jac, Xx11, Xx12, Xx13, Xx21, Xx22, Xx23, Xx31, Xx32, Xx33, 
            Jg11, Jg12, Jg13, Jg21, Jg22, Jg23, Jg31, Jg32, Jg33, nga);
}
template void gpuTemplateElemGeom3D(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
template void gpuTemplateElemGeom3D(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*,
                            float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {        
        jacg[i] = 1.0;
        nlg[i] = -1.0;  
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga, int gridDim, int blockDim)
{
    gpuTemplateFaceGeom1D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuTemplateFaceGeom1D(double*, double*, double*, int, int, int);
template void gpuTemplateFaceGeom1D(float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {
        int j = i+nga;
        jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
        nlg[i] = Jg[j]/jacg[i];
        nlg[j] = -Jg[i]/jacg[i];
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga, int gridDim, int blockDim)
{
    gpuTemplateFaceGeom2D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuTemplateFaceGeom2D(double*, double*, double*, int, int, int);
template void gpuTemplateFaceGeom2D(float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {
        int j = i+nga;
        int k = i+2*nga;
        int m = i+3*nga;
        int n = i+4*nga;
        int p = i+5*nga;
        nlg[i] = Jg[j]*Jg[p] - Jg[k]*Jg[n];
        nlg[j] = Jg[k]*Jg[m] - Jg[i]*Jg[p];
        nlg[k] = Jg[i]*Jg[n] - Jg[j]*Jg[m];
        jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
        nlg[i] = nlg[i]/jacg[i];
        nlg[j] = nlg[j]/jacg[i];
        nlg[k] = nlg[k]/jacg[i];
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuTemplateFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga, int gridDim, int blockDim)
{
    gpuTemplateFaceGeom3D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuTemplateFaceGeom3D(double*, double*, double*, int, int, int);
template void gpuTemplateFaceGeom3D(float*, float*, float*, int, int, int);

template <typename T>
__global__ void gpuTemplateApplyXx0(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int k = tid%nga; // k in [0, nga-1]                      
        fqg[tid] = udgg[tid]*Xx[k];
        tid += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuTemplateApplyXx0(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd, int gridDim, int blockDim)
{
    gpuTemplateApplyXx0<<<gridDim, blockDim>>>(fqg, udgg,  Xx, N, nga, nd);
}
template void gpuTemplateApplyXx0(double*, double*, double*, int, int, int, int, int);
template void gpuTemplateApplyXx0(float*, float*, float*, int, int, int, int, int);

template <typename T>
__global__ void gpuTemplateApplyXx1(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int k = tid%nga; // k in [0, nga-1]              
        fqg[tid] = 0.0;        
        for (int j=0; j<nd; j++)            
            fqg[tid] = fqg[tid] + udgg[tid]*Xx[k+j*nga];
        tid += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuTemplateApplyXx1(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd, int gridDim, int blockDim)
{
    gpuTemplateApplyXx1<<<gridDim, blockDim>>>(fqg, udgg,  Xx, N, nga, nd);
}
template void gpuTemplateApplyXx1(double*, double*, double*, int, int, int, int, int);
template void gpuTemplateApplyXx1(float*, float*, float*, int, int, int, int, int);


template <typename T>
__global__ void gpuTemplateApplyXx2(T *sg, T *fg, T *Xx, int N, int nga, int nd, int ncu)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int k = tid%nga; // k in [0, nga-1]              
        sg[tid] = 0.0;        
        for (int j=0; j<nd; j++)            
            sg[tid] = sg[tid] + fg[tid+nga*ncu*j]*Xx[k+j*nga];
        tid += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuTemplateApplyXx2(T *sg, T *fg,  T *Xx, int N, int nga, int nd, int ncu, int gridDim, int blockDim)
{
    gpuTemplateApplyXx2<<<gridDim, blockDim>>>(sg, fg,  Xx, N, nga, nd, ncu);
}
template void gpuTemplateApplyXx2(double*, double*, double*, int, int, int, int, int, int);
template void gpuTemplateApplyXx2(float*, float*, float*, int, int, int, int, int, int);

template <typename T>
__global__ void gpuTemplateApplyJac1(T *sg, T *jac, int N, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]              
        sg[tid] = sg[tid]*jac[n];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuTemplateApplyJac1(T *sg, T *jac, int N, int nga, int gridDim, int blockDim)
{
    gpuTemplateApplyJac1<<<gridDim, blockDim>>>(sg, jac, N, nga);
}
template void gpuTemplateApplyJac1(double*, double*, int, int, int, int);
template void gpuTemplateApplyJac1(float*, float*, int, int, int, int);

template <typename T>
__global__ void gpuTemplateApplyJac2(T *sg, T *jac, T *ug, T *fg, T *fc_u, int N, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]       
        int i = (tid-n)/nga; // i in [0, ncu-1]       
        sg[tid] = (sg[tid]+fg[tid]-ug[tid]*fc_u[i])*jac[n];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuTemplateApplyJac2(T *sg, T *jac, T *ug, T *fg, T *fc_u, int N, int nga, int gridDim, int blockDim)
{
    gpuTemplateApplyJac2<<<gridDim, blockDim>>>(sg, jac, ug, fg, fc_u, N, nga);
}
template void gpuTemplateApplyJac2(double*, double*, double*, double*, double*, int, int, int, int);
template void gpuTemplateApplyJac2(float*, float*, float*, float*, float*, int, int, int, int);

template <typename T>
__global__ void gpuTemplateApplyNormal(T *fqg, T *uhg, T *nlg, int N, int nga, int ncu, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]      
        int i = (tid-n)/nga; // i in [0, ncu-1]
        for (int j=0; j<nd; j++)            
            fqg[n+nga*i+nga*ncu*j] = uhg[n+nga*i]*nlg[n+nga*j];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuTemplateApplyNormal(T *fqg, T *uhg, T *nlg, int N, int nga, int ncu, int nd, int gridDim, int blockDim)
{
    gpuTemplateApplyNormal<<<gridDim, blockDim>>>(fqg, uhg, nlg, N, nga, ncu, nd);
}
template void gpuTemplateApplyNormal(double*, double*, double*, int, int, int, int, int, int);
template void gpuTemplateApplyNormal(float*, float*, float*, int, int, int, int,  int, int);

 template <typename T>
__global__ void gpuTemplateApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int N, int nga, int ncu, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]      
        int i = (tid-n)/nga; // i in [0, ncu-1]
        for (int j=0; j<nd; j++)            
            fqg[n+nga*i+nga*ncu*j] = uhg[n+nga*i]*nlg[n+nga*j]*jac[n];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuTemplateApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int N, int nga, int ncu, int nd, int gridDim, int blockDim)
{
    gpuTemplateApplyJacNormal<<<gridDim, blockDim>>>(fqg, uhg, nlg, jac, N, nga, ncu, nd);
}
template void gpuTemplateApplyJacNormal(double*, double*, double*, double*, int, int, int, int, int, int);
template void gpuTemplateApplyJacNormal(float*, float*, float*, float*, int, int, int, int,  int, int);


#endif

