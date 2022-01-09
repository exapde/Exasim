#ifndef __GPUARRAYPERMUTE
#define __GPUARRAYPERMUTE

template <typename T> __global__ void gpuTemplateKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = M1*M2;
    int N = N1*N2;
    int P = M*N;
    while (idx<P) {
        int i = idx%M;
        int j = (idx-i)/M;
        int ib = i%M2;
        int ia = (i-ib)/M2;
        int jb = j%N2;
        int ja = (j-jb)/N2;
        C[idx] = A[ia+M1*ja]*B[ib+M2*jb];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2)
{                
    int BLOCKDIM = 256;
    int gridDim = (M1*M2*N1*N2 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateKron<<<gridDim, BLOCKDIM>>>(C, A, B, M1, N1, M2, N2);
}

__global__ void gpuKernelIndexPermute12(int *index, int I1, int I2, int I3)
{                 
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    while (idx<N) {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        index[idx] = j+I2*i+M*k;
        idx += blockDim.x * gridDim.x;
    }
}

void gpuIndexPermute12(int *index, int I1, int I2, int I3)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelIndexPermute12<<<gridDim, BLOCKDIM>>>(index, I1, I2, I3);
}

__global__ void gpuKernelIndexPermute13(int *index, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i3+I3*i2+I3*I2*i1+N*i4;
        idx += blockDim.x * gridDim.x;
    }
}

void gpuIndexPermute13(int *index, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3*I4 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelIndexPermute13<<<gridDim, BLOCKDIM>>>(index, I1, I2, I3, I4);
}

__global__ void gpuKernelIndexPermute23(int *index, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i1+I1*i3+I1*I3*i2+N*i4;
        idx += blockDim.x * gridDim.x;
    }
}

void gpuIndexPermute23(int *index, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3*I4 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelIndexPermute23<<<gridDim, BLOCKDIM>>>(index, I1, I2, I3, I4);
}

template <typename T> __global__ void gpuTemplatePermute(T *B, T *A, int *index, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index      
    while (idx<N) {
        B[index[idx]] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuPermute(T *B, T *A, int *index, int N)
{                
    int BLOCKDIM = 256;
    int gridDim = (N + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermute<<<gridDim, BLOCKDIM>>>(B, A, index, N);
}

template <typename T> __global__ void gpuTemplatePermuteSharedMem(T *B, T *A, int *index, int N)
{
    int tx = threadIdx.x;
    int nx = blockDim.x * gridDim.x;    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index      
    __shared__ T Ashared[1025];     
    while (idx<N) {
        Ashared[index[tx]] = A[idx]; 
        __syncthreads();
        B[idx] = Ashared[tx];
        idx += nx;
    }
}

template <typename T> void gpuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM)
{                    
    int gridDim = (N + BLOCKDIM - 1) / BLOCKDIM;
    //gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermuteSharedMem<<<gridDim, BLOCKDIM>>>(B, A, index, N);
}

template <typename T> __global__ void gpuTemplatePermute12(T *B, T *A, int I1, int I2, int I3)
{            
    //int tx = threadIdx.x;                   // block thread index    
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    
    while (idx<N) {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        B[j+I2*i+M*k] = A[idx];
    }
}

template <typename T> void gpuPermute12(T *B, T *A, int I1, int I2, int I3)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermute12<<<gridDim, BLOCKDIM>>>(B, A, I1, I2, I3);
}


template <typename T> __global__ void gpuTemplatePermute13(T *B, T *A, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        B[i3+I3*i2+I3*I2*i1+N*i4] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuPermute13(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermute13<<<gridDim, BLOCKDIM>>>(B, A, I1, I2, I3, I4);
}

template <typename T> __global__ void gpuTemplatePermute23(T *B, T *A, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;              
        B[i1+I1*i3+I1*I3*i2+N*i4] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuPermute23(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermute23<<<gridDim, BLOCKDIM>>>(B, A, I1, I2, I3, I4);
}

//void gpuIndexPermute12(int*, int, int, int);
//void gpuIndexPermute13(int*, int, int, int, int);
//void gpuIndexPermute23(int*, int, int, int, int);

template void gpuKron(double*, double*, double*, int, int, int, int);
template void gpuKron(float*, float*, float*, int, int, int, int);

template void gpuPermute(double*, double*, int*, int);
template void gpuPermute(float*, float*, int*, int);
template void gpuPermuteSharedMem(double*, double*, int*, int, int);
template void gpuPermuteSharedMem(float*, float*, int*, int, int);

template void gpuPermute12(double*, double*, int, int, int);
template void gpuPermute12(float*, float*, int, int, int);
template void gpuPermute13(double*, double*, int, int, int, int);
template void gpuPermute13(float*, float*, int, int, int, int);
template void gpuPermute23(double*, double*, int, int, int, int);
template void gpuPermute23(float*, float*, int, int, int, int);

#endif
