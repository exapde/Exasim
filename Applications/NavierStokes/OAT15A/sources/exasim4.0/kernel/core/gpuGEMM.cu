#ifndef __GPUGEMM
#define __GPUGEMM

#define BLOCKDIMGEMM 256

void print3iarray(int* a, int m, int n, int p)
{
    for (int k=0; k<p; k++) {
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++)
                printf("%d   ", a[k*n*m+j*m+i]);                  
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

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

/*
    while (idx<N) {
        B[idx] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
*/
    
    //__shared__ T Ashared[256];        
    while (idx<N) {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        B[j+I2*i+M*k] = A[idx];

/*
        int l = idx%M;
        int i = l%I2;
        int j = (l-i)/I2;
        int k = (idx-l)/M;        
        B[idx] = A[j+I1*i+M*k];
*/        
        //Ashared[tx] = A[idx];
        //__syncthreads();
        //int k = (tx-l)/M;        
        //B[idx] = Ashared[j+I2*i+M*k];
        //B[idx] = A[idx];
        idx += blockDim.x * gridDim.x;
    }

}


template <typename T> void gpuPermute12(T *B, T *A, int I1, int I2, int I3)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePermute12<<<gridDim, BLOCKDIM>>>(B, A, I1, I2, I3);
}

template <int K, typename T> __global__ void gpuTemplateGEMM(T *C, T *A, T *B, int I, int J)
{            
    int tx = threadIdx.x;                   // block thread index    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    

    // load data from global memory to shared memory
    __shared__ T Ashared[BLOCKDIMGEMM];        
    if (tx<I*K)          
      Ashared[tx] = A[tx];    
            
    while (idx<J) {        
        // load data from global memory to shared memory
        __shared__ T Bshared[K*BLOCKDIMGEMM];
        
        int i, k;
        int n = K*tx;
        int p = K*idx;
        int q = I*idx;
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
       
        // thread synchronization
        __syncthreads();

        T Csub;
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+q] = Csub;                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

// C[I*J] = A[I*K] x B[K*J]
template <int K, typename T> void gpuGEMM(T *C, T *A, T *B, int I, int J)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGEMM<K><<<gridDim, BLOCKDIMGEMM>>>(C, A, B, I, J);
}

template <int N1, int N2, typename T> void gpuGEMM2DTensor(T *C, T *A1, T *A2, T *B, T *Ctmp, int M1, int M2, int K)
{
    /* B = (A1 x A2) B = A2*B*A1' = (A1*(A2*B)')'
       INPUT:  A1[M1 by N1], A2[M2 by N2], B[N1 by N2 by K] 
       OUTPUT: B[M1 by M2 by K] 
    */

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGEMM<N1>(Ctmp, A2, B, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK] 
    gpuPermute12(C, Ctmp, M2, N1, K);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGEMM<N2>(Ctmp, A1, C, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuPermute12(C, Ctmp, M1, M2, K);   
}

template <int N1, int N2, typename T> 
void gpuTensorGEMM2D(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index1, int *index2, int M1, int M2, int K)
{
    /* B = (A1 x A2) B = A2*B*A1' = (A1*(A2*B)')'
       INPUT:  A1[M1 by N1], A2[M2 by N2], B[N1 by N2 by K] 
       OUTPUT: B[M1 by M2 by K] 
    */

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGEMM<N1>(Ctmp, A2, B, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK] 
    gpuPermute(C, Ctmp, index1, M2*N1*K);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGEMM<N2>(Ctmp, A1, C, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuPermute(C, Ctmp, index2, M1*M2*K);    
}

template <int N1, int N2, typename T> 
void gpuTensorGEMM2D(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int M2, int K)
{
    /* B = (A1 x A2) B = A2*B*A1' = (A1*(A2*B)')'
       INPUT:  A1[M1 by N1], A2[M2 by N2], B[N1 by N2 by K] 
       OUTPUT: B[M1 by M2 by K] 
    */         

    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGEMM<N1>(Ctmp, A2, B, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    //int indx[BLOCKDIM];
    //cudaMemcpy(indx, index, BLOCKDIM*sizeof(int), cudaMemcpyDeviceToHost);
    //print3iarray(indx, M2, N1, BLOCKDIM/(M2*N1));

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGEMM<N2>(Ctmp, A1, C, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    if (M2*M1<17) 
        BLOCKDIM = M2*M1*32; 
    else 
        BLOCKDIM = M2*M1*16;
    if (M2*M1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <int N1, int N2, int N3, typename T> void gpuGEMM3DTensor(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int M1, int M2, int M3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGEMM<N3>(Ctmp, A3, B, M3, N2*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK] 
    gpuPermute12(C, Ctmp, M3, N2, N1*K);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGEMM<N2>(Ctmp, A2, C, M2, M3*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK] 
    gpuPermute12(C, Ctmp, M2*M3, N1, K);

    // Ctmp := A1*B -> Ctmp[M1xM2xM3xK] := A1[M1xN1]*C[N1xM2xM3xK] 
    gpuGEMM<N3>(Ctmp, A1, C, M1, M2*M3*K);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    gpuPermute12(B, Ctmp, M1, M2*M3, K);
    gpuPermute12(C, B, M2, M3, M1*K);
    //gpuPermute13(C, Ctmp, M1, M2, M3, K);
}

template <int N1, int N2, int N3, typename T> 
void gpuTensorGEMM3D(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index1, int *index2, int *index3, int M1, int M2, int M3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGEMM<N3>(Ctmp, A3, B, M3, N2*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    gpuPermute(C, Ctmp, index1, M3*N2*N1*K);    

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGEMM<N2>(Ctmp, A2, C, M2, M3*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK]     
    gpuPermute(C, Ctmp, index2, M2*M3*N1*K);    

    // Ctmp := A1*B -> Ctmp[M1xM2xM3xK] := A1[M1xN1]*C[N1xM2xM3xK] 
    gpuGEMM<N3>(Ctmp, A1, C, M1, M2*M3*K);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    gpuPermute(C, Ctmp, index3, M1*M2*M3*K);       
    //gpuPermute13(C, Ctmp, M1, M2, M3, K);
}

template <int N1, int N2, int N3, typename T> 
void gpuTensorGEMM3D(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int M2, int M3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGEMM<N3>(Ctmp, A3, B, M3, N2*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGEMM<N2>(Ctmp, A2, C, M2, M3*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK]     
    //gpuPermute12(C, Ctmp, M2*M3, N1, K);
    int N = M2*M3*N1;
    if (N<33) 
        BLOCKDIM = N*32; 
    else if (N<65)
        BLOCKDIM = N*16;
    else if (N<129)
        BLOCKDIM = N*8;
    else if (N<257)
        BLOCKDIM = N*4;
    else if (N<513)
        BLOCKDIM = N*2;
    else
        BLOCKDIM = N;       
    if (N==64) BLOCKDIM = 256;   
    if (N==216) BLOCKDIM = 864;   
    if (N==343) BLOCKDIM = 686;   
    if (N==512) BLOCKDIM = 1024;     
    gpuIndexPermute12(index, M2*M3, N1, BLOCKDIM/N);
    gpuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);

    // Ctmp := A1*B -> Ctmp[M1xM2xM3xK] := A1[M1xN1]*C[N1xM2xM3xK] 
    gpuGEMM<N3>(Ctmp, A1, C, M1, M2*M3*K);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    //gpuPermute13(C, Ctmp, M1, M2, M3, K);    
    //gpuPermute12(B, Ctmp, M1, M2*M3, K);
    //gpuPermute12(C, B, M2, M3, M1*K);
    if (N<33) 
        BLOCKDIM = N*32; 
    else if (N<65)
        BLOCKDIM = N*16;
    else if (N<129)
        BLOCKDIM = N*8;
    else if (N<257)
        BLOCKDIM = N*4;
    else if (N<513)
        BLOCKDIM = N*2;
    else
        BLOCKDIM = N;       
    if (N==64) BLOCKDIM = 256;   
    if (N==216) BLOCKDIM = 864;   
    if (N==343) BLOCKDIM = 686;   
    if (N==512) BLOCKDIM = 1024;       
    gpuIndexPermute13(index, M1, M2, M3, BLOCKDIM/N);
    gpuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);
}

template <int K, typename T> __global__ void gpuTemplateGEMM(T *C, T *A, T *B, T alpha, T beta, int I, int J)
{            
    int tx = threadIdx.x;                   // block thread index    
    //int bx = blockIdx.x;
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    
    //int P = I*K;                            // number of elements of A    

    // load data from global memory to shared memory
    __shared__ T Ashared[BLOCKDIMGEMM];        
    if (tx<I*K)          
      Ashared[tx] = A[tx];    
            
    while (idx<J) {        
        // load data from global memory to shared memory
        __shared__ T Bshared[K*BLOCKDIMGEMM];
        //__shared__ T Cshared[K*BLOCKDIMGEMM];
        
        int i, k;
        int n = K*tx;
        int p = K*idx;
        int q = I*idx;
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
            //Bshared[tx+BLOCKDIMGEMM*k] = B[tx + BLOCKDIMGEMM*k + bx*K*BLOCKDIMGEMM];         
       
        // thread synchronization
        __syncthreads();

        T Csub;
        //#pragma unroll
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            //#pragma unroll            
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+q] = alpha*Csub + beta*C[i+q];                    
            //Cshared[i+tx*I] = alpha*Csub;                        
        }

        // thread synchronization
        //__syncthreads();

        //p = bx*K*BLOCKDIMGEMM;
        //for (i=0; i<I; i++) {
        //    n = tx+i*BLOCKDIMGEMM;
        //    C[n+p] = Cshared[n] + beta*C[n+p];
        //}

        idx += blockDim.x * gridDim.x;
    }
}

// C[I*J] = A[I*K] x B[K*J]
template <int K, typename T> void gpuGEMM(T *C, T *A, T *B, T alpha, T beta, int I, int J)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGEMM<K><<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J);
}


//template<int> void gpuGEMM(double *, double *, double *, double, double, int, int);
//template<int> void gpuGEMM(float *, float *, float *, float, float, int, int);

template <int K> __global__ void gpuTemplateDGEMM(double *C, double *A, double *B, double alpha, double beta, int I, int J)
{            
    int tx = threadIdx.x;                   // block thread index    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    
    int P = I*K;                            // number of elements of A    

    // load data from global memory to shared memory
    __shared__ double Ashared[BLOCKDIMGEMM];        
    if (tx<P)          
      Ashared[tx] = A[tx];    
            
    while (idx<J) {        
        // load data from global memory to shared memory
        __shared__ double Bshared[K*BLOCKDIMGEMM];
        
        int i, k;
        int n = K*tx;
        int p = K*idx;
        int q = I*idx;
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];

        // thread synchronization
        __syncthreads();

        double Csub;
        #pragma unroll
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            #pragma unroll            
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+q] = alpha*Csub + beta*C[i+q];                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

// C[I*J] = A[I*K] x B[K*J]
template <int K> void gpuDGEMM(double *C, double *A, double *B, double alpha, double beta, int I, int J)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateDGEMM<K><<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J);
}

__global__ void gpuTemplateDGEMM2(double *C, double *A, double *B, double alpha, double beta, int I, int J, int K)
{            
    
    int tx = threadIdx.x;                   // block thread index    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    
    int P = I*K;                            // number of elements of A    

    extern __shared__ double SharedMem[];    
    double *Ashared = &SharedMem[0];
    double *Bshared = &SharedMem[P];

    // load data from global memory to shared memory
    if (tx<P)          
      Ashared[tx] = A[tx];    
            
    while (idx<J) {                
        
        int i, k;
        int n = K*tx;
        int p = K*idx;
        int q = I*idx;

        // load data from global memory to shared memory
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];

        // thread synchronization
        __syncthreads();

        double Csub;
        #pragma unroll
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            #pragma unroll            
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+q] = alpha*Csub + beta*C[i+q];                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

void gpuDGEMM2(double *C, double *A, double *B, double alpha, double beta, int I, int J, int K)
{            
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;    
    int M = I*K+K*BLOCKDIMGEMM;
    gpuTemplateDGEMM2<<<gridDim, BLOCKDIMGEMM, M>>>(C, A, B, alpha, beta, I, J, K);
}

#endif

