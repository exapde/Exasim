#ifndef __QUADRATURE
#define __QUADRATURE

#ifdef HAVE_ONETHREAD 
void opuTensorGEMM2D(dstype *C, dstype *A1, dstype *A2, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int K)
{    
    int BLOCKDIM, nn;    

    //gpuGEMM<N2>(Ctmp, A2, B, M2, N1*K);
    nn = N1*K;
    cpuGEMM(&chn, &chn, &M2, &nn, &N2, &one, A2, &M2, B, &N2, &zero, Ctmp, &M2);   
//     printArray2D(A2, M2, N2, 0);    
//     printArray2D(B, N2, N1*4, 0);   
//     printArray2D(Ctmp, M2, N1*4, 0);   
    
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    opuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));    
    opuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);
    
    //printArray2D(C, M2, N1*4, 0);      
    //printArray3D(index, M2, N1, N1*4, 0);   
    
    //gpuGEMM<N1>(Ctmp, A1, C, M1, M2*K);
    nn = M2*K;
    cpuGEMM(&chn, &chn, &M1, &nn, &N1, &one, A1, &M1, C, &N1, &zero, Ctmp, &M1);   

    if (M2*M1<17) 
        BLOCKDIM = M2*M1*32; 
    else 
        BLOCKDIM = M2*M1*16;
    if (M2*M1==64) BLOCKDIM = 256;       
    opuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    opuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);   
    
    //printArray2D(A1, M1, N1, 0);
    //printArray2D(A2, M2, N2, 0);    
    //error("here");
}

void opuTensorGEMM3D(dstype *C, dstype *A1, dstype *A2, dstype *A3, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K) 
{
    int BLOCKDIM, nn;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    //gpuGEMM<N3>(Ctmp, A3, B, M3, N2*N1*K);
    nn = N2*N1*K;
    cpuGEMM(&chn, &chn, &M3, &nn, &N3, &one, A3, &M3, B, &N3, &zero, Ctmp, &M3);   
    
    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    opuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    opuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    //gpuGEMM<N2>(Ctmp, A2, C, M2, M3*N1*K);
    nn = M3*N1*K;
    cpuGEMM(&chn, &chn, &M2, &nn, &N2, &one, A2, &M2, C, &N2, &zero, Ctmp, &M2);   

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
    opuIndexPermute12(index, M2*M3, N1, BLOCKDIM/N);
    opuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);

    // Ctmp := A1*B -> Ctmp[M1xM2xM3xK] := A1[M1xN1]*C[N1xM2xM3xK] 
    //gpuGEMM<N1>(Ctmp, A1, C, M1, M2*M3*K);
    nn = M2*M3*K;
    cpuGEMM(&chn, &chn, &M1, &nn, &N1, &one, A1, &M1, C, &N1, &zero, Ctmp, &M1);   
    
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
    opuIndexPermute13(index, M1, M2, M3, BLOCKDIM/N);
    opuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);
}

void opuTensorGEMM(dstype *C, dstype *A1, dstype *A2, dstype *A3, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K) 
{
    if (M3*N3==0)     
        opuTensorGEMM2D(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);                       
    else    
        opuTensorGEMM3D(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);                       
}
#endif       


#ifdef HAVE_OPENMP        
void cpuTensorGEMM2D(dstype *C, dstype *A1, dstype *A2, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int K)
{    
    int BLOCKDIM, nn;    

    //gpuGEMM<N2>(Ctmp, A2, B, M2, N1*K);
    nn = N1*K;
    cpuGEMM(&chn, &chn, &M2, &nn, &N2, &one, A2, &M2, B, &N2, &zero, Ctmp, &M2);   
    
    if (M2*N1<17) 
        BLOCKDIM = M2*N1*32; 
    else 
        BLOCKDIM = M2*N1*16;
    if (M2*N1==64) BLOCKDIM = 256;       
    cpuIndexPermute12(index, M2, N1, BLOCKDIM/(M2*N1));
    cpuPermuteSharedMem(C, Ctmp, index, M2*N1*K, BLOCKDIM);
    
    //gpuGEMM<N1>(Ctmp, A1, C, M1, M2*K);
    nn = M2*K;
    cpuGEMM(&chn, &chn, &M1, &nn, &N1, &one, A1, &M1, Ctmp, &N1, &zero, C, &M1);   

    if (M2*M1<17) 
        BLOCKDIM = M2*M1*32; 
    else 
        BLOCKDIM = M2*M1*16;
    if (M2*M1==64) BLOCKDIM = 256;       
    cpuIndexPermute12(index, M2, M1, BLOCKDIM/(M2*M1));
    cpuPermuteSharedMem(C, Ctmp, index, M1*M2*K, BLOCKDIM);     
}

void cpuTensorGEMM3D(dstype *C, dstype *A1, dstype *A2, dstype *A3, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K) 
{
    int BLOCKDIM, nn;

    // Ctmp := A3*B -> Ctmp[M3xN2xN1xK] := A3[M3xN3]*B[N3xN2xN1xK] 
    //gpuGEMM<N3>(Ctmp, A3, B, M3, N2*N1*K);
    nn = N2*N1*K;
    cpuGEMM(&chn, &chn, &M3, &nn, &N3, &one, A3, &M3, B, &N3, &zero, Ctmp, &M3);   
    
    // C = permute(Ctmp, [2 1 3]) -> Ctmp[N2xM3xN1*K] := C[M3xN2xN1xK]     
    if (M3*N2<17) 
        BLOCKDIM = M3*N2*32; 
    else 
        BLOCKDIM = M3*N2*16;
    if (M3*N2==64) BLOCKDIM = 256;       
    cpuIndexPermute12(index, M3, N2, BLOCKDIM/(M3*N2));
    cpuPermuteSharedMem(C, Ctmp, index, M3*N2*N1*K, BLOCKDIM);

    // Ctmp = A2*C -> Ctmp[M2xM3xN1xK] = A2[M2xN2]*C[N2xM3xN1xK] 
    //gpuGEMM<N2>(Ctmp, A2, C, M2, M3*N1*K);
    nn = M3*N1*K;
    cpuGEMM(&chn, &chn, &M2, &nn, &N2, &one, A2, &M2, B, &N2, &zero, Ctmp, &M2);   

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
    cpuIndexPermute12(index, M2*M3, N1, BLOCKDIM/N);
    cpuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);

    // Ctmp := A1*B -> Ctmp[M1xM2xM3xK] := A1[M1xN1]*C[N1xM2xM3xK] 
    //gpuGEMM<N1>(Ctmp, A1, C, M1, M2*M3*K);
    nn = M2*N3*K;
    cpuGEMM(&chn, &chn, &M1, &nn, &N1, &one, A2, &M1, B, &N1, &zero, Ctmp, &M1);   
    
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
    cpuIndexPermute13(index, M1, M2, M3, BLOCKDIM/N);
    cpuPermuteSharedMem(C, Ctmp, index, N*K, BLOCKDIM);
}

void cpuTensorGEMM(dstype *C, dstype *A1, dstype *A2, dstype *A3, dstype *B, dstype *Ctmp, 
        int *index, int M1, int N1, int M2, int N2, int M3, int N3, int K) 
{
    if (M3*N3==0)     
        cpuTensorGEMM2D(C, A1, A2, B, Ctmp, index, M1, N1, M2, N2, K);                       
    else    
        cpuTensorGEMM3D(C, A1, A2, A3, B, Ctmp, index, M1, N1, M2, N2, M3, N3, K);                       
}
#endif       

static void Node2Gauss(cublasHandle_t handle, dstype *ug, dstype *un, dstype *shapt, Int ng, Int np, Int nn, 
        dstype *ut, dstype *shapt1, dstype *shapt2, dstype *shapt3, Int *index, Int ng1, Int np1,
        Int ng2, Int np2, Int ng3, Int np3, Int nd, Int ts, Int backend)
{   
    if ((ts==1) && (nd>=2)) {     
        #ifdef HAVE_ONETHREAD         
            if (backend == 0) 
                opuTensorGEMM(ug, shapt1, shapt2, shapt3, un, ut, index, ng1, np1, ng2, np2, ng3*(nd-2), np3*(nd-2), nn);                
        #endif     
            
        #ifdef HAVE_OPENMP         
            if (backend == 1) 
                cpuTensorGEMM(ug, shapt1, shapt2, shapt3, un, ut, index, ng1, np1, ng2, np2, ng3*(nd-2), np3*(nd-2), nn);                
        #endif     
                    
        #ifdef HAVE_CUDA       
            if (backend == 2)  
                gpuTensorGEMM(ug, shapt1, shapt2, shapt3, un, ut, index, ng1, np1, ng2, np2, ng3*(nd-2), np3*(nd-2), nn);        
        #endif        
            
//         printArray3D(index, ng1, np1, 256/(ng1*np1), backend);   
//         printArray2D(shapt1, ng1, np1, backend);
//         printArray2D(shapt2, ng2, np2, backend);
//         error("here");
    }    
    else {        
        if (backend <= 1) 
            cpuGEMM(&chn, &chn, &ng, &nn, &np, &one, shapt, &ng, un, &np, &zero, ug, &ng);   

        #ifdef HAVE_CUDA          
            if (backend == 2) {
                if ((np<=48) && (ng*np<=1600))
                    gpuGEMM(ug, shapt, un, ng, nn, np);
                else
                    CHECK_CUBLAS(cublasGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, ng, nn, np, 
                        cublasOne, shapt, ng, un, np, cublasZero, ug, ng));                                
            }
        #endif                     
    }
}

static void Gauss2Node(cublasHandle_t handle, dstype *un, dstype *ug, dstype *shapg, Int ng, Int np, Int nn, 
        dstype *ut, dstype *shapg1, dstype *shapg2, dstype *shapg3, Int *index, Int ng1, Int np1, 
        Int ng2, Int np2, Int ng3, Int np3, Int nd, Int ts, Int backend)
{   
    if ((ts==1) && (nd>=2)) {               
        #ifdef HAVE_ONETHREAD    
            if (backend == 0) 
                opuTensorGEMM(un, shapg1, shapg2, shapg3, ug, ut, index, np1, ng1, np2, ng2, np3*(nd-2), ng3*(nd-2), nn);                
        #endif        
        
        #ifdef HAVE_OPENMP    
            if (backend == 1) 
                cpuTensorGEMM(un, shapg1, shapg2, shapg3, ug, ut, index, np1, ng1, np2, ng2, np3*(nd-2), ng3*(nd-2), nn);                
        #endif        
            
        #ifdef HAVE_CUDA       
            if (backend == 2)  
                gpuTensorGEMM(un, shapg1, shapg2, shapg3, ug, ut, index, np1, ng1, np2, ng2, np3*(nd-2), ng3*(nd-2), nn);                
        #endif        
    }    
    else {        
        if (backend <= 1)             
            cpuGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &zero, un, &np);  
            
        #ifdef HAVE_CUDA          
            if (backend == 2) {
                if ((ng<=48) && (ng*np<=1600))
                    gpuGEMM(un, shapg, ug, np, nn, ng);
                else
                    CHECK_CUBLAS(cublasGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_N, np, nn, ng, 
                        cublasOne, shapg, np, ug, ng, cublasZero, un, np));                    
            }
        #endif                     
    }
}
        
#endif    