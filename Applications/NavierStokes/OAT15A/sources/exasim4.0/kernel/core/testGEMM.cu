/* GEMM is a General Matrix Multiply - a subroutine in the Basic Linear Algebra Subprograms library*/

/* Includes, system */
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <cstring>
#include <sys/time.h>

/* Includes, cuda */
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "helper_cuda.h"
#include "gpuGEMM.cu"
#include "gpuArrayGemm.cu"

//using namespace std;

#define BLOCK_SIZE 16

/* ======================================================= */
/* CUDA implementation of dGEMM without using shared memory
/* ======================================================= */
__global__ void cuda_dgemm(int M, int N, int K, 
			   double alpha, 
			   const double *A, 
			   const double *B,
			   double beta, 
			   double *C) {

  int row = blockDim.y * blockIdx.y + threadIdx.y;
  int col = blockDim.x * blockIdx.x + threadIdx.x;
  
  if (row >= M || col >= N) return;
  
  double prod = 0;
  for (int k = 0; k < K; ++k){ 
    prod += A[k * M + row] * B[col * K + k];
  }
  C[col*M + row] = alpha * prod + beta * C[col*M + row]; 
}

/* ======================================================= */
/* CUDA implementation of dGEMM using shared memory
/* ======================================================= */
__global__ void cuda_dgemm_shmem(int n, 
			   double alpha, 
			   const double *B, 
			   const double *A,
			   double beta, 
			   double *C) {
  // Block index
  int block_col = blockIdx.x;
  int block_row = blockIdx.y;

  // Thread index
  int thread_col = threadIdx.x;
  int thread_row = threadIdx.y;

  //printf("row = %d col = %d  n= %d\n", block_col, block_row, n);
  //int row = blockDim.y * blockIdx.y + threadIdx.y;
  //int col = blockDim.x * blockIdx.x + threadIdx.x;
  
  int aBegin = n * blockDim.x * block_row;
  int aEnd = aBegin + n-1;
  int bBegin = blockDim.x * block_col;
  int bStep = n * blockDim.x;
  double Csub = 0;

  for (int a=aBegin, b=bBegin, istep=0;
       a <= aEnd; a+= blockDim.x, b+=bStep, ++istep){

    __shared__ double As[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE];

    if ((istep*blockDim.x+thread_col < n) && (block_row*blockDim.x+ thread_row < n))
      As[thread_row][thread_col] = A[a + n * thread_row + thread_col];
    else
      As[thread_row][thread_col] = 0;
      
    if ((block_col*blockDim.x+thread_col < n) && (istep*blockDim.x + thread_row < n))
      Bs[thread_row][thread_col] = B[b + n * thread_row + thread_col];
    else
      Bs[thread_row][thread_col] = 0;

    __syncthreads();

    // calculate the cell
    for (int k = 0; k < blockDim.x; ++k)
      Csub += As[thread_row][k] * Bs[k][thread_col];

    __syncthreads();
  }

  // Write the block sub-matrix to global memory;
  // each thread writes one element
  int c = n * blockDim.x * block_row + blockDim.x * block_col;
  if ((block_col*blockDim.x+thread_col < n) && (block_row*blockDim.x+ thread_row < n))
    C[c + n * thread_row + thread_col] = alpha * Csub + beta * C[c +n * thread_row + thread_col];

 }


/* ======================================================= */
/* Simple host implementation of a simple version of sgemm */
/* ======================================================= */
static void simple_dgemm(int M, int N, int K, double alpha, const double *A, const double *B,
                         double beta, double *C) {
  int i, j, k;
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j){
      double prod = 0;
      for (k = 0; k < K; ++k){
	prod += A[k * M + i] * B[j * K + k];
      }
      C[j * M + i] = alpha * prod + beta * C[j * M + i];
    }
  }
}

/* ======================================================= */
/* Simple host implementation of a simple version of sgemm */
/* ======================================================= */
static void simple_permute12(int M, int N, int K, double *A) 
{
  //double B[M*N*K];
  double *B = (double *)malloc(M*N*K * sizeof(double) );
  int i, j, k;
  for (k = 0; k < K; ++k) 
    for (j = 0; j < N; ++j)
      for (i = 0; i < M; ++i)
        B[i + M*j + M*N*k]  =  A[j + N*i + N*M*k];
  k = M*N*K;            
  for (i=0; i<k; i++)
    A[i] = B[i];
  free(B);
}

/* ======================= */
/* dgemm from BLAS library */
/* ======================= */
extern "C"{
extern void dgemm_(char *, char * , 
		  int *, int *, int *,
		  double *, double *, int *,
		  double *, int *,
		   double *, double *, int *); };

/* ==== */
/* Main */
/* ==== */
int main(int argc, char **argv)
{
  cublasStatus_t status;
  double *h_A, *h_B, *h_C, *h_C_blas, *h_C_simple, *h_C_0;
  double *d_A = 0; 
  double *d_B = 0;
  double *d_C = 0;
  double *d_Ct = 0;
  double alpha = 1.0f;
  double beta = 0.0f;
  int nA, nB, nC, N, M, K;
  int i;  
  cublasHandle_t handle;
  struct timeval tv1, tv2;

  M = 4; 
  K = 4;
  N = K*1024*64;
  
  nA = M * K;
  nB = K * N;
  nC = M * N;

  printf("\nRunning dgemm test for (%d by %d ) time (%d by %d) matricies.\n", M, K, K, N);
  /* Initialize CUBLAS */
  status = cublasCreate(&handle);
  
  /* Allocate host memory for the matrices */
  h_A = (double *)malloc(nA * sizeof(double) );
  h_B = (double *)malloc(nB * sizeof(double) );
  h_C = (double *)malloc(nC * sizeof(double) );
  h_C_blas = (double *)malloc(nC * sizeof(double) );
  h_C_simple = (double *)malloc(nC * sizeof(double) );
  h_C_0 = (double *)malloc(nC * sizeof(double) );

  /* Fill the matrices with test data */
  for (i = 0; i < nA; i++)
    h_A[i] = rand() / (double)RAND_MAX;

  for (i = 0; i < nB; i++)
    h_B[i] = rand() / (double)RAND_MAX;

  for (i = 0; i < nC; i++){
    h_C[i] = rand() / (double)RAND_MAX;
    h_C_blas[i] = h_C[i];
    h_C_simple[i] = h_C[i];
    h_C_0[i] = h_C[i];    
  }

  printf("\tTesting dgemm function from cuBLAS library.\n");  

  simple_dgemm(M, N, K, alpha, h_A, h_B, beta, h_C_simple);

  /* Allocate device memory for the matrices */
  cudaMalloc((void **)&d_A, nA * sizeof(d_A[0]));
  cudaMalloc((void **)&d_B, nB * sizeof(d_B[0]));
  cudaMalloc((void **)&d_C, nC * sizeof(d_C[0]));
  cudaMalloc((void **)&d_Ct, nC * sizeof(d_C[0]));

  /* Initialize the device matrices with the host matrices */
  status = cublasSetVector(nA, sizeof(h_A[0]), h_A, 1, d_A, 1);
  status = cublasSetVector(nB, sizeof(h_B[0]), h_B, 1, d_B, 1);
  status = cublasSetVector(nC, sizeof(h_C[0]), h_C, 1, d_C, 1);

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL); 
  /* Performs operation using cublas */
  status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M, N, K, &alpha, d_A, M, d_B, K, &beta, d_C, M);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  /* Read the result back */
  status = cublasGetVector(nC, sizeof(h_C[0]), d_C, 1, h_C, 1); 

  double e=0.0;
  for (int i=0; i<nC; i++)
    if (fabs(h_C[i]-h_C_simple[i])>e)
        e = fabs(h_C[i]-h_C_simple[i]);
  printf("Maximum error: %g\n", e);  

  /* ============ CUDA implementation using shared memory =============== */
  printf("\tTesting CUDA dgemm function using Shared memory.\n");  

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);  
  gpuGEMM<4>(d_C, d_A, d_B, alpha, beta, M, N);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  /* Read the result back */
  cudaMemcpy(h_C, d_C, nC*sizeof(d_C[0]), cudaMemcpyDeviceToHost);

  e=0.0;
  for (int i=0; i<nC; i++)
    if (fabs(h_C[i]-h_C_simple[i])>e)
        e = fabs(h_C[i]-h_C_simple[i]);
  printf("Maximum error: %g\n", e);  



/*
  printf("\tTesting CUDA dgemm function without using Shared memory.\n");  

  dim3 dimBlock(256, 4);
  dim3 dimGrid(N/256+1, M/4+1);

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);  
  cuda_dgemm<<<dimGrid, dimBlock>>>(M, N, K, alpha, d_A, d_B, beta, d_C);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  cudaMemcpy(h_C, d_C, nC*sizeof(d_C[0]), cudaMemcpyDeviceToHost);

  e=0.0;
  for (int i=0; i<nC; i++)
    if (fabs(h_C[i]-h_C_simple[i])>e)
        e = fabs(h_C[i]-h_C_simple[i]);
  printf("Maximum error: %g\n", e);  
*/


/*
  printf("\tTesting CUDA permute12 function.\n");  
  simple_permute12(M, K, N/K, h_C_simple);

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);  
  gpuPermute12(d_Ct, d_C, M, K, N/K);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);

  cudaMemcpy(h_C, d_Ct, nC*sizeof(d_C[0]), cudaMemcpyDeviceToHost);

  e=0.0;
  for (int i=0; i<nC; i++)
    if (fabs(h_C[i]-h_C_simple[i])>e)
        e = fabs(h_C[i]-h_C_simple[i]);
  printf("Maximum error: %g\n", e);  
*/

/*
  printf("\tTesting CUDA permute12 function.\n"); 

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);  
  gpuPermute12(d_B, d_C, M, M, N/M);
    cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t\tdone...\n");
  printf("\t\tExecution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);
*/  
 
//   /* ============ CUDA implementation using shared memory =============== */
//   printf("\tTesting CUDA dgemm function using Shared memory.\n");  
//   
//   cudaDeviceSynchronize();
//    gettimeofday(&tv1, NULL);
//   /* Kernel */
//   cuda_dgemm_shmem<<<dimGrid, dimBlock>>>(N, alpha, d_A, d_B, beta, d_C);
//   /* wait until all threads finish their job */
//   cudaDeviceSynchronize();
//   gettimeofday(&tv2, NULL);
// 
//   printf("\t\tdone...\n");
//   printf("\t\tExecution time (in millisec): %.2f\n",
// 	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
// 	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);
// 
//   /* Read the result back */
//   cudaMemcpy(h_C, d_C,n2*sizeof(d_C[0]), cudaMemcpyDeviceToHost);
// 
//   e=0.0;
//   for (int i=0; i<n2; i++)
//     if (fabs(h_C[i]-h_C_simple[i])>e)
//         e = fabs(h_C[i]-h_C_simple[i]);
//   printf("Maximum error: %g\n", e);  

  /* free cuda memory */
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
  cudaFree(d_Ct);

  /* Memory clean up */
  free(h_A);
  free(h_B);
  free(h_C);
  free(h_C_simple);
  free(h_C_blas);

  /* Shutdown */
  status = cublasDestroy(handle);

  return(0);
}
