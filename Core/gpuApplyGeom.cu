#ifndef __GPUAPPLYGEOM
#define __GPUAPPLYGEOM

template <typename T>
__global__ void  gpuTemplateApplyXx0(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int k = tid%nga; // k in [0, nga-1]                      
        fqg[tid] = udgg[tid]*Xx[k];
        tid += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuApplyXx0(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx0<<<gridDim, blockDim>>>(fqg, udgg,  Xx, N, nga, nd);
}
template void gpuApplyXx0(double*, double*, double*, int, int, int);
template void gpuApplyXx0(float*, float*, float*, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyXx1(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
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
void gpuApplyXx1(T *fqg, T *udgg,  T *Xx, int N, int nga, int nd)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx1<<<gridDim, blockDim>>>(fqg, udgg,  Xx, N, nga, nd);
}
template void gpuApplyXx1(double*, double*, double*, int, int, int);
template void gpuApplyXx1(float*, float*, float*, int, int, int);


template <typename T>
__global__ void  gpuTemplateApplyXx2(T *sg, T *fg, T *Xx, int N, int nga, int nd, int ncu)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int k = tid%nga; // k in [0, nga-1]              
        sg[tid] = 0.0;        
        for (int j=0; j<nd; j++)            
            sg[tid] = sg[tid] + fg[tid+N*j]*Xx[k+j*nga];
        tid += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuApplyXx2(T *sg, T *fg,  T *Xx, int N, int nga, int nd, int ncu)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx2<<<gridDim, blockDim>>>(sg, fg,  Xx, N, nga, nd, ncu);
}
template void gpuApplyXx2(double*, double*, double*, int, int, int, int);
template void gpuApplyXx2(float*, float*, float*, int, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyXx3(T *sg, T *ug, T *Xx, int nge, int nd, int M, int N, int P, int I, int J, int K)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {      
        int i = idx%M;       // [1, nge*ne]         
        int k = (idx-i)/M;   // [1, ncu]
        int g = i%nge;       // [1, nge]
        int e = (i-g)/nge;   // [1, ne]
        for (int m=0; m<nd; m++)
            for (int j=0; j<nd; j++)
                sg[g+nge*j+I*k+J*m+K*e] = ug[idx]*Xx[g+nge*e+M*m+P*j];
        idx += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuApplyXx3(T *sg, T *ug,  T *Xx, int nge, int nd, int ncu, int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*nd;
    int J = I*ncu;
    int K = J*nd;

    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx3<<<gridDim, blockDim>>>(sg, ug,  Xx, nge, nd, M, N, P, I, J, K);
}
template void gpuApplyXx3(double*, double*, double*, int, int, int, int);
template void gpuApplyXx3(float*, float*, float*, int, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyXx4(T *rg, T *sg, T *fg, T *Xx, T *jac, int nge, int nd, int M, int N, int P, int I, int J)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {          
        int i = idx%M;       // [1, nge*ne]         
        int k = (idx-i)/M;   // [1, ncu]
        int g = i%nge;       // [1, nge]
        int e = (i-g)/nge;   // [1, ne]
        // idx = g + nge*e + nge*ne*k
        rg[g+nge*0+I*k+J*e] = sg[idx]*jac[i];        
        for (int m=0; m<nd; m++) {
            rg[g+nge*(m+1)+I*k+J*e] = fg[idx+N*0]*Xx[g+nge*e+M*0+P*m];
            for (int j=1; j<nd; j++)
                rg[g+nge*(m+1)+I*k+J*e] += fg[idx+N*j]*Xx[g+nge*e+M*j+P*m];
        }
        idx += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuApplyXx4(T *rg, T *sg, T *fg, T *Xx, T *jac, int nge, int nd, int ncu, int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*(nd+1);
    int J = I*ncu;

    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx4<<<gridDim, blockDim>>>(rg, sg, fg, Xx, jac, nge, nd, M, N, P, I, J);
}
template void gpuApplyXx4(double*, double*, double*, double*, double*, int, int, int, int);
template void gpuApplyXx4(float*, float*, float*, float*, float*, int, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyXx5(T *rg, T *fg, T *Xx, int nge, int nd, int M, int N, int P, int I, int J)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {          
        int i = idx%M;       // [1, nge*ne]         
        int k = (idx-i)/M;   // [1, ncu]
        int g = i%nge;       // [1, nge]
        int e = (i-g)/nge;   // [1, ne]
        // idx = g + nge*e + nge*ne*k
        for (int m=0; m<nd; m++) {
            rg[g+nge*m+I*k+J*e] = fg[idx+N*0]*Xx[g+nge*e+M*0+P*m];
            for (int j=1; j<nd; j++)
                rg[g+nge*m+I*k+J*e] += fg[idx+N*j]*Xx[g+nge*e+M*j+P*m];
        }
        idx += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuApplyXx5(T *rg, T *fg, T *Xx, int nge, int nd, int ncu, int ne)
{
    int M = nge*ne;
    int N = M*ncu;
    int P = M*nd;
    int I = nge*nd;
    int J = I*ncu;    

    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyXx5<<<gridDim, blockDim>>>(rg, fg, Xx, nge, nd, M, N, P, I, J);
}
template void gpuApplyXx5(double*, double*, double*, int, int, int, int);
template void gpuApplyXx5(float*, float*, float*, int, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyJac(T *sg, T *fhg, T *jac, int nga, int ncu, int ngf, int M, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {  
        int i = idx%nga;   
        int m = (idx-i)/nga; // [1, ncu]
        int g = i%ngf;       // [1, ngf]
        int e = (i-g)/ngf;   // [1, nf]                
        sg[g+ngf*m+M*e] = fhg[idx]*jac[i];
        idx += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuApplyJac(T *sg, T *fhg, T *jac, int nga, int ncu, int ngf)
{
    int M = ngf*ncu;
    int N = nga*ncu;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJac<<<gridDim, blockDim>>>(sg, fhg, jac, nga, ncu, ngf, M, N);
}
template void gpuApplyJac(double*, double*, double*, int, int, int);
template void gpuApplyJac(float*, float*, float*, int, int, int);

template <typename T>
__global__ void  gpuTemplateApplyJac1(T *sg, T *jac, int N, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]              
        sg[tid] = sg[tid]*jac[n];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuApplyJac1(T *sg, T *jac, int N, int nga)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJac1<<<gridDim, blockDim>>>(sg, jac, N, nga);
}
template void gpuApplyJac1(double*, double*, int, int);
template void gpuApplyJac1(float*, float*, int, int);

template <typename T>
__global__ void  gpuTemplateApplyJac2(T *sg, T *jac, T *ug, T *fg, T *fc_u, int N, int nga)
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
void gpuApplyJac2(T *sg, T *jac, T *ug, T *fg, T *fc_u, int N, int nga)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJac2<<<gridDim, blockDim>>>(sg, jac, ug, fg, fc_u, N, nga);
}
template void gpuApplyJac2(double*, double*, double*, double*, double*, int, int);
template void gpuApplyJac2(float*, float*, float*, float*, float*, int, int);

template <typename T>
__global__ void  gpuTemplateApplyNormal(T *fqg, T *uhg, T *nlg, int N, int nga, int ncu, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]      
        int i = (tid-n)/nga; // i in [0, ncu-1]
        for (int j=0; j<nd; j++)            
            fqg[n+nga*i+N*j] = uhg[n+nga*i]*nlg[n+nga*j];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuApplyNormal(T *fqg, T *uhg, T *nlg, int N, int nga, int ncu, int nd)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyNormal<<<gridDim, blockDim>>>(fqg, uhg, nlg, N, nga, ncu, nd);
}
template void gpuApplyNormal(double*, double*, double*, int, int, int, int);
template void gpuApplyNormal(float*, float*, float*, int, int, int, int);

 template <typename T>
__global__ void  gpuTemplateApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int N, int nga, int ncu, int nd)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {  
        int n = tid%nga;     // n in [0, nga-1]      
        int i = (tid-n)/nga; // i in [0, ncu-1]
        for (int j=0; j<nd; j++)            
            fqg[n+nga*i+N*j] = uhg[n+nga*i]*nlg[n+nga*j]*jac[n];
        tid += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd)
{
    int N = nga*ncu;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJacNormal<<<gridDim, blockDim>>>(fqg, uhg, nlg, jac, N, nga, ncu, nd);
}
template void gpuApplyJacNormal(double*, double*, double*, double*, int, int, int);
template void gpuApplyJacNormal(float*, float*, float*, float*, int, int, int);

 template <typename T> 
__global__ void  gpuTemplateApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int N, int M, int P, int nga, int ncu, int nd, int ngf)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {  
        int i = idx%nga;   
        int m = (idx-i)/nga; // [1, ncu]
        int g = i%ngf;       // [1, ngf]
        int e = (i-g)/ngf;   // [1, nf]
        for (int j=0; j<nd; j++)
            fqg[g+ngf*m+M*j+P*e] = uhg[idx]*nlg[i+nga*j]*jac[i];
        idx += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuApplyJacNormal(T *fqg, T *uhg, T *nlg, T *jac, int nga, int ncu, int nd, int ngf)
{
    int N = nga*ncu;
    int M = ngf*ncu;
    int P = M*nd; 
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJacNormal<<<gridDim, blockDim>>>(fqg, uhg, nlg, jac, N, M, P, nga, ncu, nd, ngf);
}
template void gpuApplyJacNormal(double*, double*, double*, double*, int, int, int, int);
template void gpuApplyJacNormal(float*, float*, float*, float*, int, int, int, int);


 template <typename T> 
__global__ void  gpuTemplateShapJacDetNormal(T *fqg, T *shapgt, T *nlg, T *jac, int N, int M, int P, int nga, int npf, int nd, int ngf)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {  
        int i = idx%nga;   
        int m = (idx-i)/nga; // [1, npf]
        int g = i%ngf;       // [1, ngf]
        int e = (i-g)/ngf;   // [1, nf]
        for (int j=0; j<nd; j++)
            fqg[g+ngf*m+M*j+P*e] = shapgt[g+ngf*m]*nlg[i+nga*j]*jac[i];
        idx += blockDim.x * gridDim.x;
    }             
}

template <typename T>
void gpuShapJacDetNormal(T *fqg, T *shapgt, T *nlg, T *jac, int nga, int npf, int nd, int ngf)
{
    int N = nga*npf;
    int M = ngf*npf;
    int P = M*nd; 
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateShapJacDetNormal<<<gridDim, blockDim>>>(fqg, shapgt, nlg, jac, N, M, P, nga, npf, nd, ngf);
}
template void gpuShapJacDetNormal(double*, double*, double*, double*, int, int, int, int);
template void gpuShapJacDetNormal(float*, float*, float*, float*, int, int, int, int);

template <typename T> 
__global__  void gpuTemplateApplyFactor(T *Rfac, T *R, T *fac, int npe, int M, int N)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    while (n < N) {
        int m = n%M;       // [1, npe*ncr] 
        int k = m%npe;     // [1, npe]
        int j = (m-k)/npe; // [1, ncr]
        Rfac[n] = R[n]/fac[j];
        n += blockDim.x * gridDim.x;
    }            
}

template <typename T>
void gpuApplyFactor(T *Rfac, T *R, T *fac, int npe, int M, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyFactor<<<gridDim, blockDim>>>(Rfac, R, fac, npe, M, N);
}
template void gpuApplyFactor(double*, double*, double*, int, int, int);
template void gpuApplyFactor(float*, float*, float*, int, int, int);

template <typename T> 
__global__ void gpuTemplateApplyJac(T *Rfac, T *R, T *jac, int M, int N)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    while (n < N) {          
        int m = n%M;       // [1, npe*ncr] 
        int i = (n-m)/M;   // [1, ne] 
        Rfac[n] = R[n]*jac[i];
        n += blockDim.x * gridDim.x;
    }                        
}

template <typename T>
void gpuApplyJac(T *Rfac, T *R, T *jac, int M, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJac<<<gridDim, blockDim>>>(Rfac, R, jac, M, N);
}
template void gpuApplyJac(double*, double*, double*, int, int);
template void gpuApplyJac(float*, float*, float*, int, int);

template <typename T> 
__global__ void gpuTemplateApplyJacInv(T *Rfac, T *R, T *jac, int M, int N)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    while (n < N) {          
        int m = n%M;       // [1, npe*ncr] 
        int i = (n-m)/M;   // [1, ne] 
        Rfac[n] = R[n]/jac[i];
        n += blockDim.x * gridDim.x;
    }                        
}

template <typename T>
void gpuApplyJacInv(T *Rfac, T *R, T *jac, int M, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyJacInv<<<gridDim, blockDim>>>(Rfac, R, jac, M, N);
}
template void gpuApplyJacInv(double*, double*, double*, int, int);
template void gpuApplyJacInv(float*, float*, float*, int, int);

template <typename T> 
__global__ void gpuTemplateApplyFactorJac(T *Rfac, T *R, T *fac, T *jac, int npe, int M, int N)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    while (n < N) {          
        int m = n%M;       // [1, npe*ncr] 
        int i = (n-m)/M;   // [1, ne] 
        int k = m%npe;     // [1, npe]
        int j = (m-k)/npe; // [1, ncr]
        Rfac[n] = R[n]/(fac[j]*jac[i]);
        n += blockDim.x * gridDim.x;
    }                        
}

template <typename T>
void gpuApplyFactorJac(T *Rfac, T *R, T *fac, T *jac, int npe, int M, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyFactorJac<<<gridDim, blockDim>>>(Rfac, R, fac, jac, npe, M, N);
}
template void gpuApplyFactorJac(double*, double*, double*, double*, int, int, int);
template void gpuApplyFactorJac(float*, float*, float*, float*, int, int, int);

template <typename T> 
__global__ void gpuTemplateShapJac(T *shapjac, T *shapegt, T *jac, int nge, int M, int N)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {                                      
        int l = i%M;       // [1, nge*npe]           
        int k = (i-l)/M;   // [1, ne] 
        int n = l%nge;     // [1,nge]
        int m = (l-n)/nge; // [1, npe]
        shapjac[i] = shapegt[n+nge*m]*jac[n+nge*k];
        i += blockDim.x * gridDim.x;
    }    
}

template <typename T>
void gpuShapJac(T *shapjac, T *shapegt, T *jac, int nge, int M, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateShapJac<<<gridDim, blockDim>>>(shapjac, shapegt, jac, nge, M, N);
}
template void gpuShapJac(double*, double*, double*, int, int, int);
template void gpuShapJac(float*, float*, float*, int, int, int);

template <typename T>
__global__ void  gpuTemplateShapJacInv(T *sg, T *shapegt, T *Xx, int nge, int nd, int M, int N, int P, int I, int J, int K)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {      
        int i = idx%M;       // [1, nge*ne]         
        int k = (idx-i)/M;   // [1, npe]
        int g = i%nge;       // [1, nge]
        int e = (i-g)/nge;   // [1, ne]
        for (int m=0; m<nd; m++)
            for (int j=0; j<nd; j++)
                sg[g+nge*j+I*k+J*m+K*e] = shapegt[g+nge*k]*Xx[g+nge*e+M*m+P*j];
        idx += blockDim.x * gridDim.x;
    }    
}

template <typename T> 
void gpuShapJacInv(T *sg, T *shapegt, T *Xx, int nge, int nd, int npe, int ne)
{
    int M = nge*ne;
    int N = M*npe;
    int P = M*nd;
    int I = nge*nd;
    int J = I*npe;
    int K = J*nd;
       
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateShapJacInv<<<gridDim, blockDim>>>(sg, shapegt,  Xx, nge, nd, M, N, P, I, J, K); 
}
template void gpuShapJacInv(double*, double*, double*, int, int, int, int);
template void gpuShapJacInv(float*, float*, float*, int, int, int, int);

#endif

