#ifndef __GPUELEMFACENODE
#define __GPUELEMFACENODE

template <typename T>
__global__ void  gpuTemplateGetElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc)
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
void gpuGetElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc)
{    
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetElemNodes<<<gridDim, blockDim>>>(un, u, N, nn, np, e1, nc1, nc);
}
template void gpuGetElemNodes(double*, double*, int, int, int, int, int, int);
template void gpuGetElemNodes(float*, float*, int, int, int, int, int, int);


template <typename T>
__global__ void  gpuTemplatePutElemNodes(T *u, T *un, int N, int nn, int np, int e1, int nc1, int nc)
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
void gpuPutElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePutElemNodes<<<gridDim, blockDim>>>(un, u, N, nn, np, e1, nc1, nc);
}
template void gpuPutElemNodes(double*, double*, int, int, int, int, int, int);
template void gpuPutElemNodes(float*, float*, int, int, int, int, int, int);

template <typename T>
__global__ void  gpuTemplateGetFaceNodes0(T *uh, T *udg, int *facecon, int N, int ndf, 
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
void gpuGetFaceNodes0(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetFaceNodes0<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuGetFaceNodes0(double*, double*, int*, int, int, int, int, int, int);
template void gpuGetFaceNodes0(float*, float*, int*, int, int, int, int, int, int);


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
void gpuGetFaceNodes1(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetFaceNodes1<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuGetFaceNodes1(double*, double*, int*, int, int, int, int, int, int);
template void gpuGetFaceNodes1(float*, float*, int*, int, int, int, int, int, int);


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
void gpuGetFaceNodes2(T *uh, T *udg, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetFaceNodes2<<<gridDim, blockDim>>>(uh, udg, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuGetFaceNodes2(double*, double*, int*, int, int, int, int, int, int);
template void gpuGetFaceNodes2(float*, float*, int*, int, int, int, int, int, int);

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
void gpuPutFaceNodes0(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePutFaceNodes0<<<gridDim, blockDim>>>(udg, uh, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuPutFaceNodes0(double*, double*, int*, int, int, int, int, int, int);
template void gpuPutFaceNodes0(float*, float*, int*, int, int, int, int, int, int);


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
void gpuPutFaceNodes1(T *udg, T *uh, int *facecon, int N, int ndf, 
                   int npf, int npe, int nc, int f1)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePutFaceNodes1<<<gridDim, blockDim>>>(udg, uh, facecon, N, ndf, npf, npe, nc, f1);
}
template void gpuPutFaceNodes1(double*, double*, int*, int, int, int, int, int, int);
template void gpuPutFaceNodes1(float*, float*, int*, int, int, int, int, int, int);

template <typename T>
__global__ void  gpuTemplatePutFaceNodesV0(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1, int *rowe2f2, 
    int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int M, int N, int K, int I)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;    
    while (idx < N) {
        int j = idx%M;              //[1, npe*nc]
        int k = (idx-j)/M+e1;       //[1, ne]      
        int l = j%npe;              //[1, npe]
        int m = (j-l)/npe;          //[1, nc] 
        int q, p, s;

        int i = ent2ind1[l+npe*k];
        int e = (i > 0) ? i : 0;
        int n = rowe2f1[i+1] - rowe2f1[e];
        for (j=0; j<n; j++) {
            q = cole2f1[rowe2f1[i]+j];
            p = q%npf;              // [1, npf]
            s = (q-p)/npf;          // [1, nf]    
            udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
        }            

        i = ent2ind2[l+npe*k];
        e = (i > 0) ? i : 0;
        n = rowe2f2[i+1] - rowe2f2[e];
        for (j=0; j<n; j++) {
            q = cole2f2[rowe2f2[i]+j];
            p = q%npf;              // [1, npf]
            s = (q-p)/npf;          // [1, nf]      
            udg[I+idx] = udg[I+idx] + uh[p+npf*m+K*s]; 
        }                                
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void  gpuTemplatePutFaceNodesV1(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1, int *rowe2f2, 
    int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int M, int N, int K, int I)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;    
    while (idx < N) {
        int j = idx%M;              //[1, npe*nc]
        int k = (idx-j)/M+e1;       //[1, ne]      
        int l = j%npe;              //[1, npe]
        int m = (j-l)/npe;          //[1, nc] 
        int q, p, s;

        int i = ent2ind1[l+npe*k];
        int e = (i > 0) ? i : 0;
        int n = rowe2f1[i+1] - rowe2f1[e];
        for (j=0; j<n; j++) {
            q = cole2f1[rowe2f1[i]+j];
            p = q%npf;              // [1, npf]
            s = (q-p)/npf;          // [1, nf]    
            udg[I+idx] = udg[I+idx] - uh[p+npf*m+K*s]; 
        }            
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts)
{
    int ne = e2-e1;
    int K = npf*nc;
    int M = npe*nc;
    int N = M*ne;        
    int I = M*e1;

    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    if (opts==0)
        gpuTemplatePutFaceNodesV0<<<gridDim, blockDim>>>(udg, uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, 
            cole2f2, ent2ind2, npf, npe, nc, e1, e2, M, N, K, I);
    else
        gpuTemplatePutFaceNodesV1<<<gridDim, blockDim>>>(udg, uh, rowe2f1, cole2f1, ent2ind1, rowe2f2, 
            cole2f2, ent2ind2, npf, npe, nc, e1, e2, M, N, K, I);
}
template void gpuPutFaceNodes(double*, double*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void gpuPutFaceNodes(float*, float*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

#endif

