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

void print1darray(double* a, int m)
{    
    for (int i=0; i<m; i++)
        printf("%g   ", a[i]);          
    printf("\n");
}

void print2darray(double* a, int m, int n)
{
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            printf("%g   ", a[j*m+i]);              
        printf("\n");
    }
    printf("\n");
}

void print3darray(double* a, int m, int n, int p)
{
    for (int k=0; k<p; k++) {
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++)
                printf("%g   ", a[k*n*m+j*m+i]);                  
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

void printArray2D(double* a, int m, int n)
{
    int N = m*n;
    double *b = (double*) malloc (sizeof (double)*N);
    cudaMemcpy(b, a, N*sizeof(double), cudaMemcpyDeviceToHost);    
    print2darray(b, m, n);
    free(b);
}

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
  double *h_A, *h_B, *h_C, *h_D, *h_E, *h_F, *h_A1, *h_A2, *h_A3;
  double *d_A, *d_B, *d_C, *d_D, *d_E, *d_F, *d_A1, *d_A2, *d_A3;
  double *h_Cts, *d_Cts, *d_Ctm, *h_Fts, *d_Fts, *d_Ftm;
  int i, K, M1, N1, M2, N2, M3, N3, S1, S2, S3, SA, SB, SC, SD, SE, SF, SI;
  int *index1, *index2;    
  cublasHandle_t handle;
  struct timeval tv1, tv2;
  status = cublasCreate(&handle);
  double alpha = 1.0f;
  double beta = 0.0f;
 
  K = 1024*16;
  M1 = 8; M2 = 8; M3 = 8;
  N1 = 8; N2 = 8; N3 = 8;
  S1 = M1*N1;
  S2 = M2*N2;
  S3 = M3*N3;
  SA = S1*S2;
  SB = N1*N2*K;  
  SC = M1*M2*K;
  SD = S1*S2*S3;  
  SE = N1*N2*N3*K;  
  SF = M1*M2*M3*K;  
  SI = M2*N1*K;

  h_A1 = (double *)malloc(S1 * sizeof(double) );  
  h_A2 = (double *)malloc(S2 * sizeof(double) );  
  h_A3 = (double *)malloc(S3 * sizeof(double) );  
  h_A = (double *)malloc(SA * sizeof(double) );  
  h_B = (double *)malloc(SB * sizeof(double) );  
  h_C = (double *)malloc(SC * sizeof(double) );  
  h_D = (double *)malloc(SD * sizeof(double) );  
  h_E = (double *)malloc(SE * sizeof(double) );  
  h_F = (double *)malloc(SF * sizeof(double) );  
  h_Cts = (double *)malloc(SC * sizeof(double) );  
  h_Fts = (double *)malloc(SF * sizeof(double) );

  for (i = 0; i < S1; i++)
    h_A1[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < S2; i++)
    h_A2[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < S3; i++)
    h_A3[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SA; i++)
    h_A[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SB; i++)
    h_B[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SC; i++)
    h_C[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SD; i++)
    h_D[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SE; i++)
    h_E[i] = rand() / (double)RAND_MAX;
  for (i = 0; i < SF; i++)
    h_F[i] = rand() / (double)RAND_MAX;

  cudaMalloc((void **)&d_A1, S1 * sizeof(double));
  cudaMalloc((void **)&d_A2, S2 * sizeof(double));
  cudaMalloc((void **)&d_A3, S3 * sizeof(double));
  cudaMalloc((void **)&d_A, SA * sizeof(double));
  cudaMalloc((void **)&d_B, SB * sizeof(double));
  cudaMalloc((void **)&d_C, SC * sizeof(double));
  cudaMalloc((void **)&d_D, SD * sizeof(double));
  cudaMalloc((void **)&d_E, SE * sizeof(double));
  cudaMalloc((void **)&d_F, SF * sizeof(double));
  cudaMalloc((void **)&d_Cts, SC * sizeof(double));
  cudaMalloc((void **)&d_Ctm, SC * sizeof(double));
  cudaMalloc((void **)&d_Fts, SF * sizeof(double));
  cudaMalloc((void **)&d_Ftm, SF * sizeof(double));  
  cudaMalloc((void **)&index1, SI * sizeof(int));
  cudaMalloc((void **)&index2, SC * sizeof(int));

  gpuIndexPermute12(index1, M2, N1, K);
  gpuIndexPermute12(index2, M1, M2, K);

  status = cublasSetVector(S1, sizeof(h_A1[0]), h_A1, 1, d_A1, 1);
  status = cublasSetVector(S2, sizeof(h_A2[0]), h_A2, 1, d_A2, 1);
  status = cublasSetVector(S3, sizeof(h_A3[0]), h_A3, 1, d_A3, 1);
  status = cublasSetVector(SA, sizeof(h_A[0]), h_A, 1, d_A, 1);
  status = cublasSetVector(SB, sizeof(h_B[0]), h_B, 1, d_B, 1);
  status = cublasSetVector(SC, sizeof(h_C[0]), h_C, 1, d_C, 1);
  status = cublasSetVector(SD, sizeof(h_D[0]), h_D, 1, d_D, 1);
  status = cublasSetVector(SE, sizeof(h_E[0]), h_E, 1, d_E, 1);
  status = cublasSetVector(SF, sizeof(h_F[0]), h_F, 1, d_F, 1);

  gpuKron(d_A, d_A1, d_A2, M1, N1, M2, N2);  
  gpuKron(d_D, d_A, d_A3, M1*M2, N1*N2, M3, N3);

/*
  print2darray(h_A1, M1, N1);
  printArray2D(d_A1, M1, N1);
  print2darray(h_A2, M2, N2);
  printArray2D(d_A2, M2, N2);
  print2darray(h_A3, M3, N3);
  printArray2D(d_A3, M3, N3);
  printArray2D(d_A, M1*M2, N1*N2);
  printArray2D(d_D, M1*M2*M3, N1*N2*N3);
*/
  
  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);   
  status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M1*M2, K, N1*N2, &alpha, d_A, M1*M2, d_B, N1*N2, &beta, d_C, M1*M2);
  //status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M2, N1*K, N2, &alpha, d_A2, M2, d_B, N2, &beta, d_Ctm, M2);
  //printArray2D(d_Ctm, M2, N1); 
  //gpuPermute12(d_C, d_Ctm, M2, N1, K);
  //printArray2D(d_C, N1, M2);
  //status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M1, M2*K, N1, &alpha, d_A1, M1, d_C, N1, &beta, d_Ctm, M1);
  //printArray2D(d_Ctm, M1, M2);
  //gpuPermute12(d_C, d_Ctm, M1, M2, K);   
  //printArray2D(d_C, M2, M1);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t 2D cublasDgemm execution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);  

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);   
  //gpuGEMM2DTensor<8,8>(d_Cts, d_A1, d_A2, d_B, d_Ctm, M1, M2, K);
  //gpuTensorGEMM2D<8,8>(d_Cts, d_A1, d_A2, d_B, d_Ctm, index1, index2, M1, M2, K);
  gpuTensorGEMM2D<8,8>(d_Cts, d_A1, d_A2, d_B, d_Ctm, index1, M1, M2, K);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t 2D tensor execution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);  

/*
  printArray2D(d_A1, M1, N1);
  printArray2D(d_A2, M2, N2);
  printArray2D(d_A, M1*M2, N1*N2);
  printArray2D(d_B, N1, N2);

  printArray2D(d_C, M1, M2);
  printArray2D(d_Cts, M1, M2);
*/

  status = cublasGetVector(SC, sizeof(h_C[0]), d_C, 1, h_C, 1); 
  status = cublasGetVector(SC, sizeof(h_C[0]), d_Cts, 1, h_Cts, 1);
  double e=0.0;
  for (int i=0; i<M1*M2*K; i++)
    if (fabs(h_C[i]-h_Cts[i])>e)
        e = fabs(h_C[i]-h_Cts[i]);
  printf("Maximum error: %g\n", e);  
  
  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);   
  status = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M1*M2*M3, K, N1*N2*N3, &alpha, d_D, M1*M2*M3, d_E, N1*N2*N3, &beta, d_F, M1*M2*M3);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t 3D cublasDgemm execution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);  

  cudaDeviceSynchronize();
  gettimeofday(&tv1, NULL);   
  //gpuGEMM3DTensor<8,8,8>(d_Fts, d_A1, d_A2, d_A3, d_E, d_Ftm, M1, M2, M3, K);
  //gpuTensorGEMM3D(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, int M1, int M2, int M3, int K);
  gpuTensorGEMM3D<8,8,8>(d_Fts, d_A1, d_A2, d_A3, d_E, d_Ftm, index1, M1, M2, M3, K);
  cudaDeviceSynchronize();
  gettimeofday(&tv2, NULL);

  printf("\t 3D tensor execution time (in millisec): %.2f\n",
	 (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
	 (double)(tv2.tv_sec -tv1.tv_sec )*1000);  

  status = cublasGetVector(SF, sizeof(h_F[0]), d_F, 1, h_F, 1); 
  status = cublasGetVector(SF, sizeof(h_F[0]), d_Fts, 1, h_Fts, 1);
  e=0.0;
  for (int i=0; i<M1*M2*M3*K; i++)
    if (fabs(h_F[i]-h_Fts[i])>e)
        e = fabs(h_F[i]-h_Fts[i]);
  printf("Maximum error: %g\n", e);    

  cudaFree(d_A1);
  cudaFree(d_A2);
  cudaFree(d_A3);
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
  cudaFree(d_D);
  cudaFree(d_E);
  cudaFree(d_F);
  cudaFree(d_Cts);
  cudaFree(d_Ctm);
  
  free(h_A1);
  free(h_A2);
  free(h_A3);
  free(h_A);
  free(h_B);
  free(h_C);
  free(h_D);
  free(h_E);
  free(h_F);
  free(h_Cts);

  /* Shutdown */
  status = cublasDestroy(handle);

  return(0);
}
