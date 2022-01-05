#ifndef __GPUARRAYGEMM
#define __GPUARRAYGEMM

#define BLOCKDIMGEMM 256

template <typename T> 
__global__ void gpuTemplateGemmV0(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K, int N)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {        
        int i = idx%I;
        int j = (idx-i)/I;                
        int m = K*j;        
        T prod = 0.0;
        for (int k=0; k<K; k++)
            prod += A[i+I*k]*B[k+m];
        C[idx] = alpha*prod + beta*C[idx];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuGemmV0(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int N = I*J;        
    int gridDim = (N + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV0<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K, N);
}

template <typename T> 
__global__ void gpuTemplateGemmV1(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K, int M, int N)
{        
    // static shared memory
    __shared__ T Ashared[BLOCKDIMGEMM];

    // load data from global memory to shared memory
    int idx = threadIdx.x;
    while (idx<M) {          
      Ashared[idx] = alpha*A[idx];    
      idx = idx + blockDim.x;
    }

    // thread synchronization
    __syncthreads();

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {        
        int i = idx%I;
        int j = (idx-i)/I;                
        int m = K*j;
        T prod = 0.0;        
        for (int k=0; k<K; k++)
            prod += Ashared[i+I*k]*B[k+m];
        C[idx] = prod + beta*C[idx];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuGemmV1(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int N = I*J;        
    int gridDim = (N + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV1<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K, I*K, N);
}

template <typename T> 
__global__ void gpuTemplateGemmV2(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K, int M, int N)
{        
    // static shared memory
    __shared__ T Ashared[1280];

    // load data from global memory to shared memory
    int idx = threadIdx.x;
    while (idx<M) {          
      Ashared[idx] = alpha*A[idx];    
      idx = idx + blockDim.x;
    }

    // thread synchronization
    __syncthreads();

    idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {        
        int i = idx%I;
        int j = (idx-i)/I;                
        int m = K*j;
        T prod = 0.0;        
        for (int k=0; k<K; k++)
            prod += Ashared[i+I*k]*B[k+m];
        C[idx] = prod + beta*C[idx];        
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuGemmV2(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int N = I*J;        
    int gridDim = (N + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV2<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K, I*K, N);
}

template<typename T> __global__ void gpuTemplateGemmV3(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx<J) {                        
        int m = K*idx;        
        int p = I*idx;
        T Csub;
        for (int i=0; i<I; i++) {               
            Csub = beta*C[i+p];
            for (int k=0; k<K; k++)
                Csub += alpha*A[i+I*k]*B[k+m];
            C[i+p] = Csub;                    
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGemmV3(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV3<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K);
}

template<typename T> __global__ void gpuTemplateGemmV4(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K, int M)
{                
    // load data from global memory to shared memory
    __shared__ T Ashared[BLOCKDIMGEMM];            
    int idx = threadIdx.x;
    while (idx<M) {          
      Ashared[idx] = alpha*A[idx];    
      idx = idx + blockDim.x;
    }

    // thread synchronization
    __syncthreads();
            
    idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index    
    while (idx<J) {                        
        int m = K*idx;        
        int p = I*idx;
        T Csub;
        for (int i=0; i<I; i++) {               
            Csub = beta*C[i+p];
            for (int k=0; k<K; k++)
                Csub += Ashared[i+I*k]*B[k+m];
            C[i+p] = Csub;                    
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGemmV4(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV4<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K, I*K);
}

template<typename T> __global__ void gpuTemplateGemmV5(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K, int M)
{                
    __shared__ T Ashared[1280];

    int idx = threadIdx.x;
    while (idx<M) {          
      Ashared[idx] = alpha*A[idx];    
      idx = idx + blockDim.x;
    }

    // thread synchronization
    __syncthreads();
            
    idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index    
    while (idx<J) {                        
        int m = K*idx;        
        int p = I*idx;
        T Csub;
        for (int i=0; i<I; i++) {               
            Csub = beta*C[i+p];
            for (int k=0; k<K; k++)
                Csub += Ashared[i+I*k]*B[k+m];
            C[i+p] = Csub;                    
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGemmV5(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV5<<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, K, I*K);
}

template <int K, typename T> __global__ void gpuTemplateGemmV6(T *C, T *A, T *B, T alpha, T beta, int I, int J, int M)
{            
    int tx = threadIdx.x;                   // block thread index        

    // load data from global memory to shared memory
    __shared__ T Ashared[BLOCKDIMGEMM];            
    int idx = tx;
    while (idx<M) {          
      Ashared[idx] = alpha*A[idx];    
      idx = idx + blockDim.x;
    }

    idx = tx + blockIdx.x * blockDim.x; // global thread index             
    while (idx<J) {        
        // load data from global memory to shared memory
        __shared__ T Bshared[K*BLOCKDIMGEMM];
        
        int i, k;
        int n = K*tx;
        int p = K*idx;        
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
       
        // thread synchronization
        __syncthreads();

        p = I*idx;
        T Csub;
        for (i=0; i<I; i++) {               
            Csub = beta*C[i+p];
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+p] = Csub;                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

template <int K, typename T> void gpuGemmV6(T *C, T *A, T *B, T alpha, T beta, int I, int J)
{                
    int gridDim = (J + BLOCKDIMGEMM - 1) / BLOCKDIMGEMM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGemmV6<K><<<gridDim, BLOCKDIMGEMM>>>(C, A, B, alpha, beta, I, J, I*K);
}

template <typename T> 
__global__ void gpuTemplateGemmV7(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{            
  int row = blockDim.y * blockIdx.y + threadIdx.y;
  int col = blockDim.x * blockIdx.x + threadIdx.x;
    
  while (col < J) {
    T prod = 0;
    for (int k = 0; k < K; ++k) 
        prod += A[k * I + row] * B[col * K + k];  
    C[col*I + row] = alpha * prod + beta * C[col*I + row]; 
    
    col += blockDim.x * gridDim.x;
  }
}

template <typename T> void gpuGemmV7(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int BDIMX;
    if (I<3)
        BDIMX = 128;
    else if (I<6)
        BDIMX = 64;
    else if (I<16)
        BDIMX = 32;
    else if (I<32)
        BDIMX = 16;
    else if (I<64)
        BDIMX = 8;
    else if (I<128)
        BDIMX = 4;
    else
        BDIMX = 2;

    dim3 block (BDIMX, I);
    int GDIMX = (J + block.x - 1) / block.x;
    GDIMX = (GDIMX>1024)? 1024 : GDIMX;
    dim3 grid (GDIMX, 1);
    gpuTemplateGemmV7<<<grid, block>>>(C, A, B, alpha, beta, I, J, K);
}

template <typename T> 
__global__ void gpuTemplateGemmV8(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{            
  int row = blockDim.y * blockIdx.y + threadIdx.y;
  int col = blockDim.x * blockIdx.x + threadIdx.x;
    
  __shared__ T Ashared[BLOCKDIMGEMM];                
  if (col<K)           
      Ashared[col * I + row] = A[col * I + row];    

  __syncthreads();

  while (col < J) {
    T prod = 0;
    for (int k = 0; k < K; ++k) 
        prod += Ashared[k * I + row] * B[col * K + k];  
    C[col*I + row] = alpha * prod + beta * C[col*I + row]; 
    
    col += blockDim.x * gridDim.x;
  }
}

template <typename T> void gpuGemmV8(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int BDIMX;
    if (I<3)
        BDIMX = 128;
    else if (I<6)
        BDIMX = 64;
    else if (I<16)
        BDIMX = 32;
    else if (I<32)
        BDIMX = 16;
    else if (I<64)
        BDIMX = 8;
    else if (I<128)
        BDIMX = 4;
    else
        BDIMX = 2;

    dim3 block (BDIMX, I);
    int GDIMX = (J + block.x - 1) / block.x;
    GDIMX = (GDIMX>1024)? 1024 : GDIMX;
    dim3 grid (GDIMX, 1);
    gpuTemplateGemmV8<<<grid, block>>>(C, A, B, alpha, beta, I, J, K);
}

template <typename T> 
__global__ void gpuTemplateGemmV9(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{            
  int row = blockDim.y * blockIdx.y + threadIdx.y;
  int col = blockDim.x * blockIdx.x + threadIdx.x;
    
  __shared__ T Ashared[1280];                
  while (col<K) {          
      Ashared[col * I + row] = A[col * I + row];    
      col += blockDim.x * gridDim.x;
  }
  __syncthreads();

  col = blockDim.x * blockIdx.x + threadIdx.x;
  while (col < J) {
    T prod = 0;
    for (int k = 0; k < K; ++k) 
        prod += Ashared[k * I + row] * B[col * K + k];  
    C[col*I + row] = alpha * prod + beta * C[col*I + row]; 
    
    col += blockDim.x * gridDim.x;
  }
}

template <typename T> void gpuGemmV9(T *C, T *A, T *B, T alpha, T beta, int I, int J, int K)
{        
    int BDIMX;
    if (I<3)
        BDIMX = 128;
    else if (I<6)
        BDIMX = 64;
    else if (I<16)
        BDIMX = 32;
    else if (I<32)
        BDIMX = 16;
    else if (I<64)
        BDIMX = 8;
    else if (I<128)
        BDIMX = 4;
    else
        BDIMX = 2;

    dim3 block (BDIMX, I);
    int GDIMX = (J + block.x - 1) / block.x;
    GDIMX = (GDIMX>1024)? 1024 : GDIMX;
    dim3 grid (GDIMX, 1);
    gpuTemplateGemmV9<<<grid, block>>>(C, A, B, alpha, beta, I, J, K);
}


template <typename T> 
void gpuTensorGemm2DV1(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int N1, int M2, int N2, int K)
{
    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGemmV1(Ctmp, A2, B, (T) 1.0, (T) 0.0, M2, N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGemmV1(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*K, N1);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <typename T> 
void gpuTensorGemm2DV4(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int N1, int M2, int N2, int K)
{
    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGemmV4(Ctmp, A2, B, (T) 1.0, (T) 0.0, M2, N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGemmV4(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*K, N1);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <typename T> 
void gpuTensorGemm2DV8(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int N1, int M2, int N2, int K)
{
    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGemmV8(Ctmp, A2, B, (T) 1.0, (T) 0.0, M2, N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGemmV8(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*K, N1);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <int N1, int N2, typename T> 
void gpuTensorGemm2DV6(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int M2, int K)
{
    /* B = (A1 x A2) B = A2*B*A1' = (A1*(A2*B)')'
       INPUT:  A1[M1 by N1], A2[M2 by N2], B[N1 by N2 by K] 
       OUTPUT: B[M1 by M2 by K] 
    */         

    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGemmV6<N2>(Ctmp, A2, B, (T) 1.0, (T) 0.0, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;      
 
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGemmV6<N1>(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <typename T> 
void gpuTensorGemm3DV1(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGemmV1(Ctmp, A3, B, (T) 1.0, (T) 0.0, M3, N2*N1*K, N3);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGemmV1(Ctmp, A2, C, (T) 1.0, (T) 0.0, M2, M3*N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK]     
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
    gpuGemmV1(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*M3*K, N1);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    N = M2*M3*M1;
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

template <typename T> 
void gpuTensorGemm3DV4(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGemmV4(Ctmp, A3, B, (T) 1.0, (T) 0.0, M3, N2*N1*K,  N3);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGemmV4(Ctmp, A2, C, (T) 1.0, (T) 0.0, M2, M3*N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK]     
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
    gpuGemmV4(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*M3*K, N1);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    N = M2*M3*M1;
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

template <typename T> 
void gpuTensorGemm3DV8(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGemmV8(Ctmp, A3, B, (T) 1.0, (T) 0.0, M3, N2*N1*K,  N3);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGemmV8(Ctmp, A2, C, (T) 1.0, (T) 0.0, M2, M3*N1*K, N2);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xM3xK] := Ctmp[M2xM3xN1xK]     
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
    gpuGemmV8(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*M3*K, N1);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
    N = M2*M3*M1;
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

template <int N1, int N2, int N3, typename T> 
void gpuTensorGemm3DV6(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int M2, int M3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGemmV6<N3>(Ctmp, A3, B, (T) 1.0, (T) 0.0, M3, N2*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGemmV6<N2>(Ctmp, A2, C, (T) 1.0, (T) 0.0, M2, M3*N1*K);

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
    gpuGemmV6<N1>(Ctmp, A1, C, (T) 1.0, (T) 0.0, M1, M2*M3*K);

    // C = permute(Ctmp, [3 2 1 4]) -> C[M3xM2xM1xK] := Ctmp[M1xM2xM3xK] 
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

template <typename T> 
void gpuTensorGemmV1(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    if (M3*N3==0)     
        gpuTensorGemm2DV1(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);                
    else 
        gpuTensorGemm3DV1(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);                
}

template <typename T> 
void gpuTensorGemmV4(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    if (M3*N3==0)     
        gpuTensorGemm2DV4(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);                
    else 
        gpuTensorGemm3DV4(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);                
}

template <typename T> 
void gpuTensorGemmV8(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    if (M3*N3==0)     
        gpuTensorGemm2DV8(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);                
    else 
        gpuTensorGemm3DV8(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);                
}

template <typename T> 
void gpuTensorGemmV6(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    if (M3*N3==0) 
    {
        if ((N1==3) && (N2==3))        
            gpuTensorGemm2DV6<3,3>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==4) && (N2==4))        
            gpuTensorGemm2DV6<4,4>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==5) && (N2==5))        
            gpuTensorGemm2DV6<5,5>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==6) && (N2==6))        
            gpuTensorGemm2DV6<6,6>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==7) && (N2==7))        
            gpuTensorGemm2DV6<7,7>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==8) && (N2==8))        
            gpuTensorGemm2DV6<8,8>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else {            
        }           
    }        
    else
    {
        if ((N1==3) && (N2==3) && (N3==3))        
            gpuTensorGemm3DV6<3,3,3>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==4) && (N2==4) && (N3==4))        
            gpuTensorGemm3DV6<4,4,4>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==5) && (N2==5) && (N3==5))        
            gpuTensorGemm3DV6<5,5,5>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==6) && (N2==6) && (N3==6))        
            gpuTensorGemm3DV6<6,6,6>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==7) && (N2==7) && (N3==7))        
            gpuTensorGemm3DV6<7,7,7>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==8) && (N2==8) && (N3==8))        
            gpuTensorGemm3DV6<8,8,8>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else {            
        }            
    }
}

template void gpuGemmV0(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV0(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV1(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV1(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV2(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV2(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV3(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV3(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV4(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV4(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV5(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV5(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV7(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV7(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV8(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV8(float*, float*, float*, float, float, int, int, int);
template void gpuGemmV9(double*, double*, double*, double, double, int, int, int);
template void gpuGemmV9(float*, float*, float*, float, float, int, int, int);
template void gpuTensorGemmV1(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV1(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV4(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV4(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV8(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV8(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV6(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGemmV6(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);

#endif

