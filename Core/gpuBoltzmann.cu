#ifndef __GPUBOLTZMANN
#define __GPUBOLTZMANN

template <typename T> 
__global__ void gpuTemplateGetxg(T* xb, T* xdg, int* fb, int ngf, int nfe, int ne, int M, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int m = idx%M;     // [1, nbf*ngf]
        int i = (idx-m)/M; // [1, ncx]
        int j = m%ngf;     // [1, ngf]
        int k = (m-j)/ngf; // [1, nbf]
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        xb[idx] = xdg[j+ngf*f+ngf*nfe*e+ngf*nfe*ne*i];        
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpugetxg(T* xb, T* xdg, int* fb, int nbf, int ngf, int nfe, int ne, int ncx)
{        
    int M = nbf*ngf;
    int N = nbf*ngf*ncx;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetxg<<<gridDim, blockDim>>>(xb, xdg, fb, ngf, nfe, ne, M, N);
}
template void gpugetxg(double*, double*, int*, int, int, int, int, int);
template void gpugetxg(float*, float*, int*, int, int, int, int, int);

template <typename T> 
__global__ void gpuTemplateGetug(T* ug, T* udg, int* fb, int* perm, int npf, int npe, int ne, int M, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int m = idx%M;     // [1, nbf*npf]
        int i = (idx-m)/M; // [1, ncu]
        int j = m%npf;     // [1, npf]
        int k = (m-j)/npf; // [1, nbf]
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        ug[idx] = udg[perm[j+npf*f]+npe*e+npe*ne*i];         
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpugetug(T* ug, T* udg, int* fb, int* perm, int nbf, int npf, int npe, int ne, int ncu)
{        
    int M = nbf*npf;
    int N = nbf*npf*ncu;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetug<<<gridDim, blockDim>>>(ug, udg, fb, perm, npf, npe, ne, M, N);
}
template void gpugetug(double*, double*, int*, int*, int, int, int, int, int);
template void gpugetug(float*, float*, int*, int*, int, int, int, int, int);

template <typename T> 
__global__ void gpuTemplateGetfg(T* fg, T* fdg, int* fb, int* perm, int npf, int npe, int ne, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int j = idx%npf;     // [1, npf]
        int k = (idx-j)/npf; // [1, nbf]        
        int e = fb[0+4*k];
        int f = fb[1+4*k];
        int n = fb[2+4*k];
        fg[idx] = fdg[perm[j+npf*f]+npe*e+npe*ne*n];
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpugetfg(T* fg, T* fdg, int* fb, int* perm, int nbf, int npf, int npe, int ne)
{        
    int N = nbf*npf;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetfg<<<gridDim, blockDim>>>(fg, fdg, fb, perm, npf, npe, ne, N);
}
template void gpugetfg(double*, double*, int*, int*, int, int, int, int);
template void gpugetfg(float*, float*, int*, int*, int, int, int, int);

template <typename T> 
__global__ void gpuTemplatePutfg(T* fdg, T* fg, T scalar, int* fc, int M, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < M) {
        int K = fc[(N+1)+idx+1]-fc[(N+1)+idx]; 
        for (int j=0; j<K; j++)
            fdg[fc[idx]] += scalar*fg[fc[2*(N+1)+fc[(N+1)+idx]+j]];
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpuputfg(T* fdg, T* fg, T scalar, int* fc, int M, int N)
{            
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePutfg<<<gridDim, blockDim>>>(fdg, fg, scalar, fc, M, N);
}
template void gpuputfg(double*, double*, double, int*, int, int);
template void gpuputfg(float*, float*, float, int*, int, int);


template <typename T> 
__global__ void gpuTemplateGetc(T* c, T* pv, int* fb, int nd, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int j = idx%nd;     // [1, nd]
        int k = (idx-j)/nd; // [1, nbf]        
        int n = fb[2+4*k];
        c[j+nd*k] = pv[j+nd*n]; 
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpugetc(T* c, T* pv, int* fb, int nbf, int nd)
{            
    int N = nbf*nd;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetc<<<gridDim, blockDim>>>(c, pv, fb, nd, N);
}
template void gpugetc(double*, double*, int*, int, int);
template void gpugetc(float*, float*, int*, int, int);


template <typename T> 
__global__ void gpuTemplateGetnjc(T *njc, T* DetJf, T* nlf, T *c, int ngf, int nd, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int j = idx%ngf;     // [1, ngf]
        int i = (idx-j)/ngf; // [1, nbf]      
        njc[idx] = 0.0;
        for (int k=0; k<nd; k++) 
            njc[idx] += DetJf[idx]*nlf[j+ngf*i+N*k]*c[k+nd*i];        
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpugetnjc(T *njc, T* DetJf, T* nlf, T *c, int ngf, int nbf, int nd)
{            
    int N = nbf*ngf;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetnjc<<<gridDim, blockDim>>>(njc, DetJf, nlf, c, ngf, nd, N);
}
template void gpugetnjc(double*, double*, double*, double*, int, int, int);
template void gpugetnjc(float*, float*, float*, float*, int, int, int);


template <typename T> 
__global__ void gpuTemplateApplyInverse(T* F, T* R, T *f, T* Aeinv, T* Beinv, int* t2t,
        int* nvsweep, int* evsweep1, int* evsweep2, int* linflowi,  int* ninflowi,  
        int* finflowi, int npe, int nfe, int ne, int s, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int e = evsweep1[nvsweep[s]+idx];
        int n = evsweep2[nvsweep[s]+idx];
        int j, p;
        for (j=0; j<npe; j++)
            f[npe*ne*n+npe*e+j] = R[npe*ne*n+npe*e+j];

        int k = e+(ne+1)*n;            
        int n1 = linflowi[n];                      
        int n2 = npe*npe*(n1+ninflowi[k]); 
        int nfaces = ninflowi[k+1]-ninflowi[k];
        for (int m=0; m<nfaces; m++) {
            int l = finflowi[n1+ninflowi[k]+m];
            int e1 = t2t[l+nfe*e];
            for (j=0; j<npe; j++)
                for (p=0; p<npe; p++)
                    f[npe*ne*n+npe*e+j] = f[npe*ne*n+npe*e+j]-Beinv[n2+npe*npe*m+npe*p+j]*F[p+npe*e1+npe*ne*n];                
        }
        
        for (j=0; j<npe; j++) {
            F[npe*ne*n+npe*e+j] = 0.0;
            for (p=0; p<npe; p++)
                F[npe*ne*n+npe*e+j] += Aeinv[npe*npe*ne*n+npe*npe*e+npe*p+j]*f[npe*ne*n+npe*e+p];            
        }        
        idx += blockDim.x * gridDim.x;       
    }        
}
template <typename T> void gpuapplyinverse(T* F, T* R, T *f, T* Aeinv, T* Beinv, int* t2t,
        int* nvsweep, int* evsweep1, int* evsweep2, int* linflowi, int* ninflowi,  
        int* finflowi, int npe, int nfe, int ne, int N, int s)
{            
    //int N = nvsweep[s+1]-nvsweep[s];
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateApplyInverse<<<gridDim, blockDim>>>(F, R, f, Aeinv, Beinv, t2t, nvsweep, evsweep1,
            evsweep2, linflowi, ninflowi, finflowi, npe, nfe, ne, s, N);
}
template void gpuapplyinverse(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void gpuapplyinverse(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

#endif

