#ifndef __PBLAS_H__
#define __PBLAS_H__

static void cpuNode2Gauss(dstype *ug, dstype *un, dstype *shapt, Int ng, Int np, Int nn)
{    
    //dstype alpha = 1.0, beta = 0.0;
    //char chn = 'N';
    //DGEMM(&chn, &chn, &ng, &nn, &np, &alpha, shapt, &ng, un, &np, &beta, ug, &ng);      
#ifdef USE_FLOAT        
    SGEMM(&chn, &chn, &ng, &nn, &np, &one, shapt, &ng, un, &np, &zero, ug, &ng);   
#else        
    DGEMM(&chn, &chn, &ng, &nn, &np, &one, shapt, &ng, un, &np, &zero, ug, &ng);   
#endif    
}

static void cpuGauss2Node(dstype *un, dstype *ug, dstype *shapg, Int ng, Int np, Int nn)
{    
    //dstype alpha = 1.0, beta = 0.0;
    //char chn = 'N';
    //DGEMM(&chn, &chn, &np, &nn, &ng, &alpha, shapg, &np, ug, &ng, &beta, un, &np);    
#ifdef USE_FLOAT    
    SGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &zero, un, &np);    
#else    
    DGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &zero, un, &np);    
#endif    
}

static void cpuGauss2Node1(dstype *un, dstype *ug, dstype *shapg, Int ng, Int np, Int nn)
{    
    //dstype alpha = 1.0, beta = 0.0;
    //char chn = 'N';
    //DGEMM(&chn, &chn, &np, &nn, &ng, &alpha, shapg, &np, ug, &ng, &beta, un, &np);    
#ifdef USE_FLOAT        
    SGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &one, un, &np);    
#else        
    DGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &one, un, &np);    
#endif    
}


static void Node2Gauss(cublasHandle_t handle, dstype *ug, dstype *un, dstype *shapt, Int ng, Int np, Int nn, Int backend)
{            
#ifdef USE_FLOAT        
    if (backend <= 1) 
        SGEMM(&chn, &chn, &ng, &nn, &np, &one, shapt, &ng, un, &np, &zero, ug, &ng);   
#else        
    if (backend <= 1) 
        DGEMM(&chn, &chn, &ng, &nn, &np, &one, shapt, &ng, un, &np, &zero, ug, &ng);   
#endif    

#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, ng, nn, np, 
            cublasOne, shapt, ng, un, np, cublasZero, ug, ng);    
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, ng, nn, np, 
            cublasOne, shapt, ng, un, np, cublasZero, ug, ng);
#endif        
#endif             
}

static void Gauss2Node(cublasHandle_t handle, dstype *un, dstype *ug, dstype *shapg, Int ng, Int np, Int nn, Int backend)
{            
#ifdef USE_FLOAT        
    if (backend <= 1) 
        SGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &zero, un, &np);    
#else        
    if (backend <= 1) 
        DGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &zero, un, &np);    
#endif    

#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, np, nn, ng, 
            cublasOne, shapg, np, ug, ng, cublasZero, un, np);    
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, np, nn, ng, 
            cublasOne, shapg, np, ug, ng, cublasZero, un, np);
#endif        
#endif             
}

static void Gauss2Node1(cublasHandle_t handle, dstype *un, dstype *ug, dstype *shapg, Int ng, Int np, Int nn, Int backend)
{            
#ifdef USE_FLOAT        
    if (backend <= 1) 
        SGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &one, un, &np);    
#else        
    if (backend <= 1) 
        DGEMM(&chn, &chn, &np, &nn, &ng, &one, shapg, &np, ug, &ng, &one, un, &np);    
#endif    

#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, np, nn, ng, 
            cublasOne, shapg, np, ug, ng, cublasOne, un, np);    
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, np, nn, ng, 
            cublasOne, shapg, np, ug, ng, cublasOne, un, np);
#endif        
#endif             
}

#ifdef HAVE_CUDA       
static void gpuComputeInverse(cublasHandle_t handle, dstype* A, dstype *C, Int n, Int batchSize)
{
    Int *ipiv, *info;
    cudaMalloc(&ipiv, n * batchSize * sizeof(Int));
    cudaMalloc(&info,  batchSize * sizeof(Int));     
    
    //dstype *C;
    //CHECK(cudaMalloc(&C,  n*n*batchSize * sizeof(dstype)));        
    
    dstype **Ap_h = (dstype **)malloc(batchSize*sizeof(dstype *));
    dstype **Ap_d;
    cudaMalloc(&Ap_d,batchSize*sizeof(dstype *));
    Ap_h[0] = A;
    for (Int i = 1; i < batchSize; i++)
      Ap_h[i] = Ap_h[i-1]+(n*n);
    cudaMemcpy(Ap_d,Ap_h,batchSize*sizeof(dstype *),cudaMemcpyHostToDevice);
    
    dstype **Cp_h = (dstype **)malloc(batchSize*sizeof(dstype *));
    dstype **Cp_d;
    cudaMalloc(&Cp_d,batchSize*sizeof(dstype *));
    Cp_h[0] = C;
    for (Int i = 1; i < batchSize; i++)
      Cp_h[i] = Cp_h[i-1] + (n*n);
    cudaMemcpy(Cp_d,Cp_h,batchSize*sizeof(dstype *),cudaMemcpyHostToDevice);
        
#ifdef USE_FLOAT        
    cublasSgetrfBatched(handle,n,Ap_d,n,ipiv,info,batchSize);    
    cublasSgetriBatched(handle,n,(const dstype **)Ap_d,n,ipiv,Cp_d,n,info,batchSize);
#else            
    cublasDgetrfBatched(handle,n,Ap_d,n,ipiv,info,batchSize);    
    cublasDgetriBatched(handle,n,(const dstype **)Ap_d,n,ipiv,Cp_d,n,info,batchSize);    
#endif            
    
    // copy C to A
    gpuArrayCopy(A, C, n*n*batchSize);
            
    cudaFree(Ap_d); free(Ap_h);
    cudaFree(Cp_d); free(Cp_h);
    cudaFree(ipiv); cudaFree(info); 
}
#endif     

static void cpuComputeInverse(dstype* A, dstype* work, Int* ipiv, Int n)
{
    Int lwork = n*n;
    Int info;
#ifdef USE_FLOAT           
    SGETRF(&n,&n,A,&n,ipiv,&info);
    SGETRI(&n,A,&n,ipiv,work,&lwork,&info);    
#else            
    DGETRF(&n,&n,A,&n,ipiv,&info);
    DGETRI(&n,A,&n,ipiv,work,&lwork,&info);
#endif        
}
   
static void Inverse(cublasHandle_t handle, dstype* A, dstype *C, Int *ipiv, Int n, Int batchSize, Int backend)
{    
#ifdef HAVE_CUDA        
    if (backend == 2)   
        gpuComputeInverse(handle, A, C, n, batchSize);
#else       
    if (backend <= 1) {        
        for (int k=0; k<batchSize; k++) 
            cpuComputeInverse(&A[n*n*k], C, ipiv, n);                
    }    
#endif              
}

static void PDOT(cublasHandle_t handle, Int m, dstype* x, Int incx, dstype* y, Int incy, 
        dstype *global_dot, dstype *local_dot, Int backend) 
{           
#ifdef USE_FLOAT    
    if (backend <= 1) 
        *local_dot = SDOT(&m, x, &incx, y, &incy);
#else   
    if (backend <= 1) 
        *local_dot = DDOT(&m, x, &incx, y, &incy);
#endif        
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSdot(handle, m, x, incx, y, incy, local_dot);    
#else            
    if (backend == 2)  
        cublasDdot(handle, m, x, incx, y, incy, local_dot);    
#endif        
   // cudaDeviceSynchronize();
#endif             
    
#ifdef HAVE_MPI        
#ifdef USE_FLOAT        
    MPI_Allreduce(local_dot, global_dot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else        
    MPI_Allreduce(local_dot, global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif         
#else    
    //ArrayCopy(global_dot, local_dot, 1, backend);
    *global_dot = *local_dot;
#endif    
}

static dstype PNORM(cublasHandle_t handle, Int m, dstype* x, 
        dstype* global_dot, dstype* local_dot, Int backend) 
{              
    PDOT(handle, m, x, inc1, x, inc1, global_dot, local_dot, backend); 
    return sqrt(*global_dot);    
}

static void PDOT(cublasHandle_t handle, Int m, dstype* x, Int incx, dstype* y, Int incy, 
        dstype *global_dot, Int backend) 
{           
    dstype local_dot=zero;
    INIT_TIMING;
    START_TIMING;

#ifdef USE_FLOAT    
    if (backend <= 1) 
        local_dot = SDOT(&m, x, &incx, y, &incy);
#else   
    if (backend <= 1) 
        local_dot = DDOT(&m, x, &incx, y, &incy);
#endif        
    
#ifdef HAVE_CUDA  
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSdot(handle, m, x, incx, y, incy, &local_dot);    
#else            
    if (backend == 2)  
        cublasDdot(handle, m, x, incx, y, incy, &local_dot);    
#endif
    //cudaDeviceSynchronize();        
#endif             
    END_TIMING_DISC(94);
    
    START_TIMING;
#ifdef HAVE_MPI        
#ifdef USE_FLOAT        
    MPI_Allreduce(&local_dot, global_dot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else        
    MPI_Allreduce(&local_dot, global_dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif         
#else    
    //ArrayCopy(global_dot, local_dot, 1, backend);
    *global_dot = local_dot;
#endif    
    END_TIMING_DISC(95);
}

static dstype PNORM(cublasHandle_t handle, Int m, dstype* x, Int backend) 
{            
    dstype nrm;
    PDOT(handle, m, x, inc1, x, inc1, &nrm, backend); 
    return sqrt(nrm);    
}

static void DOT(cublasHandle_t handle, Int m, dstype* x, Int incx, dstype* y, Int incy, dstype *dot, Int backend) 
{               
#ifdef USE_FLOAT    
    if (backend <= 1) 
        *dot = SDOT(&m, x, &incx, y, &incy);
#else   
    if (backend <= 1) 
        *dot = DDOT(&m, x, &incx, y, &incy);
#endif        
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSdot(handle, m, x, incx, y, incy, dot);    
#else            
    if (backend == 2)  
        cublasDdot(handle, m, x, incx, y, incy, dot);    
#endif        
#endif                 
}

static dstype NORM(cublasHandle_t handle, Int m, dstype* x, Int backend) 
{            
    dstype nrm;
    DOT(handle, m, x, inc1, x, inc1, &nrm, backend); 
    return sqrt(nrm);    
}

static void PGEMNV(cublasHandle_t handle, Int m, Int n, dstype* alpha, dstype* A, Int lda, 
        dstype* x, Int incx, dstype* beta, dstype* y, Int incy, Int backend) 
{
    /* y = alpha*A * x + beta y */
#ifdef USE_FLOAT     
    if (backend <= 1) 
        SGEMV(&chn, &m, &n, alpha, A, &lda, x, &incx, beta, y, &incy);
#else   
    if (backend <= 1) 
        DGEMV(&chn, &m, &n, alpha, A, &lda, x, &incx, beta, y, &incy);
#endif
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemv(handle, CUBLAS_OP_N, m, n, alpha, A, lda, x, incx,
                             beta, y, incy);        
#else            
    if (backend == 2)  
        cublasDgemv(handle, CUBLAS_OP_N, m, n, alpha, A, lda, x, incx,
                             beta, y, incy);
#endif        
#endif                     
}

static void PGEMTV(cublasHandle_t handle, Int m, Int n, dstype *alpha, dstype* A, Int lda, 
        dstype* x, Int incx, dstype *beta, dstype* y, Int incy, dstype* ylocal, Int backend) 
{
    /* y = alpha*A^T * x + beta y */
#ifdef USE_FLOAT     
    if (backend <= 1) 
        SGEMV(&cht, &m, &n, alpha, A, &lda, x, &incx, beta, ylocal, &incy);
#else   
    if (backend <= 1) 
        DGEMV(&cht, &m, &n, alpha, A, &lda, x, &incx, beta, ylocal, &incy);
#endif
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemv(handle, CUBLAS_OP_T, m, n, alpha, A, lda, x, incx,
                             beta, ylocal, incy);        
#else            
    if (backend == 2)  
        cublasDgemv(handle, CUBLAS_OP_T, m, n, alpha, A, lda, x, incx,
                             beta, ylocal, incy);
#endif        
#endif             
    
#ifdef  HAVE_MPI          
#ifdef USE_FLOAT         
    MPI_Allreduce(ylocal, y, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else            
    MPI_Allreduce(ylocal, y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif    
#else
    ArrayCopy(y, ylocal, n, backend);
#endif    
}

static void PGEMTV2(cublasHandle_t handle, Int m, Int n, dstype *alpha, dstype* A, Int lda, 
        dstype* x, Int incx, dstype *beta, dstype* y, Int incy, dstype* ylocal, Int backend) 
{
    /* y = alpha*A^T * x + beta y */
#ifdef USE_FLOAT     
    if (backend <= 1)
        SGEMM(&chn, &chn, &incx, &n, &m, alpha, x, &incx, A, &lda, beta, ylocal, &incy);    
        //SGEMV(&cht, &m, &n, alpha, A, &lda, x, &incx, beta, ylocal, &incy);
#else   
    if (backend <= 1) 
        DGEMM(&chn, &chn, &incx, &n, &m, alpha, x, &incx, A, &lda, beta, ylocal, &incy);  
        //DGEMV(&cht, &m, &n, alpha, A, &lda, x, &incx, beta, ylocal, &incy);
#endif
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, incx, n, m, 
            alpha, x, incx, A, lda, beta, ylocal, incy);                    
        //CHECK_CUBLAS(cublasSgemv(handle, CUBLAS_OP_T, m, n, alpha, A, lda, x, incx,
        //                     beta, ylocal, incy));        
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, incx, n, m, 
            alpha, x, incx, A, lda, beta, ylocal, incy);                    
        //CHECK_CUBLAS(cublasDgemv(handle, CUBLAS_OP_T, m, n, alpha, A, lda, x, incx,
        //                     beta, ylocal, incy));
#endif        
#endif             
    
#ifdef  HAVE_MPI          
#ifdef USE_FLOAT         
    MPI_Allreduce(ylocal, y, n, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else            
    MPI_Allreduce(ylocal, y, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif    
#else
    ArrayCopy(y, ylocal, n, backend);
#endif    
}

static void PGEMNM(cublasHandle_t handle, Int m, Int n, Int k, dstype *alpha, dstype* A, Int lda, 
        dstype* B, Int ldb, dstype *beta, dstype* C, Int ldc, Int backend) 
{
    /* C = alpha*A * B + beta C */
#ifdef USE_FLOAT     
    if (backend <= 1) 
        SGEMM(&chn, &chn, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, C, &ldc);
#else   
    if (backend <= 1) 
        DGEMM(&chn, &chn, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, C, &ldc);
#endif
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, 
            alpha, A, lda, B, ldb, beta, C, ldc);                    
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, 
            alpha, A, lda, B, ldb, beta, C, ldc);
#endif        
#endif                     
}

static void PGEMTM(cublasHandle_t handle, Int m, Int n, Int k, dstype *alpha, dstype* A, Int lda, 
        dstype* B, Int ldb, dstype *beta, dstype* C, Int ldc, dstype* Clocal, Int backend) 
{
    /* C = alpha*A^T * B + beta C */
#ifdef USE_FLOAT     
    if (backend <= 1) 
        SGEMM(&cht, &chn, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, Clocal, &ldc);
#else   
    if (backend <= 1) 
        DGEMM(&cht, &chn, &m, &n, &k, alpha, A, &lda, B, &ldb, beta, Clocal, &ldc);
#endif
    
#ifdef HAVE_CUDA          
#ifdef USE_FLOAT  
    if (backend == 2)     
        cublasSgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, 
            alpha, A, lda, B, ldb, beta, Clocal, ldc);                    
#else            
    if (backend == 2)  
        cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k, 
            alpha, A, lda, B, ldb, beta, Clocal, ldc);
#endif        
#endif                     
    
    Int p = m*n;
    
#ifdef  HAVE_MPI          
#ifdef USE_FLOAT         
    MPI_Allreduce(Clocal, C, p, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else            
    MPI_Allreduce(Clocal, C, p, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif    
#else
    ArrayCopy(C, Clocal, p, backend);
#endif        
}

// static void PGEMNMStridedBached(cublasHandle_t handle, Int m, Int n, Int k, dstype *alpha, dstype* A, Int lda, 
//         dstype* B, Int ldb, dstype *beta, dstype* C, Int ldc, Int batchCount, Int backend) 
// {
//     /* C = alpha*A * B + beta C */
// //         gemmStridedBatched(cublasHandle_t handle, 
// //                       cublasOperation_t transA, cublasOperation_t transB,
// //                       int M, int N, int K, 
// //                       const T* alpha,
// //                       const T* A, int ldA, int strideA, 
// //                       const T* B, int ldB, int strideB, 
// //                       const T* beta,
// //                       T* C, int ldC, int strideC,
// //                       int batchCount)               
// #ifdef HAVE_CUDA          
// #ifdef USE_FLOAT  
//     if (backend == 2)     
//         CHECK_CUBLAS(cublasSgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, 
//             alpha, A, lda, m*k, B, ldb, k*n, beta, C, ldc, m*n, batchCount));                    
// #else            
//     if (backend == 2)  
//         CHECK_CUBLAS(cublasDgemmStridedBatched(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, 
//             alpha, A, lda, m*k, B, ldb, k*n, beta, C, ldc, m*n, batchCount));
// #endif        
// #endif                     
// }

#endif  

