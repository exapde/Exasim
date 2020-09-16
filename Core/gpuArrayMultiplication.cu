#ifndef __GPUARRAYMULTIPLICATION
#define __GPUARRAYMULTIPLICATION

#define BLOCKDIMGEMM 256

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
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
       
        // thread synchronization
        __syncthreads();

        p = I*idx;
        T Csub;
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+p] = Csub;                    
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

template <typename T> __global__ void gpuTemplateSmallSizeGEMM(T *C, T *A, T *B, int I, int J, int K)
{            
    int tx = threadIdx.x;                   // block thread index    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    
    int i, k, n, p;

    // load data from global memory to shared memory
    __shared__ T Ashared[128];        
    if (tx<I*K)          
      Ashared[tx] = A[tx];    
            
    while (idx<J) {        
        // load data from global memory to shared memory
        __shared__ T Bshared[1024];
                
        n = K*tx;
        p = K*idx;        
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
       
        // thread synchronization
        __syncthreads();

        p = I*idx;
        T Csub;
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+p] = Csub;                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

template<typename T> __global__ void gpuTemplateBigSizeGEMM(T *C, T *A, T *B, int I, int J, int K)
{            
    int tx = threadIdx.x;                   // block thread index    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index    
    int i, k;
    i = tx;
    
    __shared__ T Ashared[1600];            
    while (i<I*K) {          
      Ashared[i] = A[i];    
      i = i + blockDim.x;
    }

    while (idx<J) {                
        __shared__ T Bshared[2048];                
        int n = K*tx;
        int p = K*idx;        
        for (k=0; k<K; k++)  
            Bshared[k+n] = B[k+p];
       
        // thread synchronization
        __syncthreads();

        p = I*idx;
        T Csub;
        for (i=0; i<I; i++) {               
            Csub = 0.0;
            for (k=0; k<K; k++)
                Csub += Ashared[i+I*k]*Bshared[k+n];
            C[i+p] = Csub;                    
        }

        idx += blockDim.x * gridDim.x;
    }
}

// C[I*J] = A[I*K] x B[K*J]
template <typename T> void gpuGEMM(T *C, T *A, T *B, int I, int J, int K)
{                    
    if (K<9) {
        int BLOCKDIM = 256;
        if (K>4)
           BLOCKDIM = 128;
        int gridDim = (J + BLOCKDIM - 1) / BLOCKDIM;
        gridDim = (gridDim>1024)? 1024 : gridDim;
        gpuTemplateSmallSizeGEMM<<<gridDim, BLOCKDIM>>>(C, A, B, I, J, K);
    }
    else {
        int BLOCKDIM = 64;
        if (K>24)
            BLOCKDIM = 32;
        int gridDim = (J + BLOCKDIM - 1) / BLOCKDIM;
        gridDim = (gridDim>1024)? 1024 : gridDim;
        gpuTemplateBigSizeGEMM<<<gridDim, BLOCKDIM>>>(C, A, B, I, J, K);
    }
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
    gpuGEMM<N2>(Ctmp, A2, B, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGEMM<N1>(Ctmp, A1, C, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    if (M2*M1<17) 
        BLOCKDIM = M2*M1*32; 
    else 
        BLOCKDIM = M2*M1*16;
    if (M2*M1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

template <typename T> 
void gpuTensorGEMM2D(T *C, T *A1, T *A2, T *B, T *Ctmp, int *index, int M1, int N1, int M2, int N2, int K)
{
    /* B = (A1 x A2) B = A2*B*A1' = (A1*(A2*B)')'
       INPUT:  A1[M1 by N1], A2[M2 by N2], B[N1 by N2 by K] 
       OUTPUT: B[M1 by M2 by K] 
    */         

    int BLOCKDIM;    

    // Ctmp := A2*B -> Ctmp[M2xN1xK] := A2[M2xN2]*B[N2xN1xK] 
    gpuGEMM(Ctmp, A2, B, M2, N1*K, N2);
    //gpuGEMM<4>(Ctmp, A2, B, M2, N1*K);

    // C = permute(Ctmp, [2 1 3]) -> C[N1xM2xK] := Ctmp[M2xN1xK]            
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    gpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM1xK] = A1[M1xN1]*C[N1xM2xK] 
    gpuGEMM(Ctmp, A1, C, M1, M2*K, N1);
    //gpuGEMM<4>(Ctmp, A1, C, M1, M2*K);

    // C = permute(Ctmp, [2 1 3]) -> C[M2xM1xK] := Ctmp[M1xM2xK] 
    if (M2*M1<17) 
        BLOCKDIM = M2*M1*32; 
    else 
        BLOCKDIM = M2*M1*16;
    if (M2*M1==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    gpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}


template <typename T> 
void gpuTensorGEMM3D(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    /* B = (A1 x A2 x A3) B 
       INPUT:  A1[M1 by N1], A2[M2 by N2], A3[M3 by N3], B[N1 by N2 by N3 by K] 
       OUTPUT: B[M1 by M2 by M3 by K] 
    */
    int BLOCKDIM;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    gpuGEMM(Ctmp, A3, B, M3, N2*N1*K,  N3);
    //gpuGEMM<5>(Ctmp, A3, B, M3, N2*N1*K);

    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    gpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    gpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    gpuGEMM(Ctmp, A2, C, M2, M3*N1*K, N2);
    //gpuGEMM<5>(Ctmp, A2, C, M2, N3*N1*K);

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
    gpuGEMM(Ctmp, A1, C, M1, M2*M3*K, N1);
    //gpuGEMM<5>(Ctmp, A1, C, M1, M2*M1*K);

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
    gpuGEMM<N1>(Ctmp, A1, C, M1, M2*M3*K);

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
void gpuTensorGEMM(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K)
{
    if (M3*N3==0) 
    {
        gpuTensorGEMM2D(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);            
/*
        if ((N1==3) && (N2==3))        
            gpuTensorGEMM2D<3,3>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==4) && (N2==4))        
            gpuTensorGEMM2D<4,4>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==5) && (N2==5))        
            gpuTensorGEMM2D<5,5>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==6) && (N2==6))        
            gpuTensorGEMM2D<6,6>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==7) && (N2==7))        
            gpuTensorGEMM2D<7,7>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else if ((N1==8) && (N2==8))        
            gpuTensorGEMM2D<8,8>(C, A1, A2, B, Ctmp, index, M1, M2, K);            
        else {            
        }           
*/ 
    }        
    else
    {
        gpuTensorGEMM3D(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);            
/*
        if ((N1==3) && (N2==3) && (N3==3))        
            gpuTensorGEMM3D<3,3,3>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==4) && (N2==4) && (N3==4))        
            gpuTensorGEMM3D<4,4,4>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==5) && (N2==5) && (N3==5))        
            gpuTensorGEMM3D<5,5,5>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==6) && (N2==6) && (N3==6))        
            gpuTensorGEMM3D<6,6,6>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==7) && (N2==7) && (N3==7))        
            gpuTensorGEMM3D<7,7,7>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else if ((N1==8) && (N2==8) && (N3==8))        
            gpuTensorGEMM3D<8,8,8>(C, A1, A2, A3, B, Ctmp, index, M1, M2, M3, K);            
        else {            
        }            
*/
    }
}

template void gpuGEMM(double*, double*, double*, int, int, int);
template void gpuGEMM(float*, float*, float*, int, int, int);
template void gpuTensorGEMM2D(double*, double*, double*, double*, double*, int*, int, int, int, int, int);
template void gpuTensorGEMM2D(float*, float*, float*, float*, float*, int*, int, int, int, int, int);
template void gpuTensorGEMM3D(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGEMM3D(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);

template<int> void gpuGEMM(double*, double*, double*, int, int);
template<int> void gpuGEMM(float*, float*, float*, int, int);
template<int, int> void gpuTensorGEMM2D(double*, double*, double*, double*, double*, int*, int, int, int);
template<int, int> void gpuTensorGEMM2D(float*, float*, float*, float*, float*, int*, int, int, int);
template<int, int, int> void gpuTensorGEMM3D(double*, double*, double*, double*, double*, double*, int*, int, int, int, int);
template<int, int, int> void gpuTensorGEMM3D(float*, float*, float*, float*, float*, float*, int*, int, int, int, int);

template void gpuTensorGEMM(double*, double*, double*, double*, double*, double*, int*, int, int, int, int, int, int, int);
template void gpuTensorGEMM(float*, float*, float*, float*, float*, float*, int*, int, int, int, int, int, int, int);

#endif

