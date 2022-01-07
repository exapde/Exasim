#ifndef __GPUAPPLYFHAT
#define __GPUAPPLYFHAT

template <typename T>   
__global__ void gpuTemplateAverageFlux(T *fg, int N)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < N) {
		fg[tid+N] = 0.5*(fg[tid] + fg[tid+N]);        
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAverageFlux(T *fg, int N)
{	
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
	gpuTemplateAverageFlux<<<gridDim, blockDim>>>(fg, N);
}

template void gpuAverageFlux(double *, int);
template void gpuAverageFlux(float *, int);


template <typename T>   
__global__ void gpuTemplateAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < M) {
        int i = tid%numPoints;                
		fg[tid] = fg[N + 0*M + tid] * nl[i + 0 * numPoints];   
        for (int l = 1; l < nd; l++)
            fg[tid] += fg[N + l*M + tid] * nl[i + l * numPoints];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAverageFluxDotNormal(T *fg, T *nl, int N, int M, int numPoints, int nd)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateAverageFluxDotNormal<<<gridDim, blockDim>>>(fg, nl, N, M, numPoints, nd);
}
template void gpuAverageFluxDotNormal(double *, double *, int, int, int, int);
template void gpuAverageFluxDotNormal(float *, float *, int, int, int, int);


template <typename T>   
__global__ void gpuTemplateAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < M) {
		fg[tid] += tau[0] * (ug1[tid] - ug2[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAddStabilization1(T *fg, T *ug1, T *ug2, T *tau, int M)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateAddStabilization1<<<gridDim, blockDim>>>(fg, ug1, ug2, tau, M);
}

template void gpuAddStabilization1(double *, double *, double *, double *, int);
template void gpuAddStabilization1(float *, float *, float *, float *, int);


template <typename T>   
__global__ void gpuTemplateAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < M) {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
		fg[tid] += tau[j] * (ug1[tid] - ug2[tid]);
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAddStabilization2(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateAddStabilization2<<<gridDim, blockDim>>>(fg, ug1, ug2, tau, M, numPoints);
}

template void gpuAddStabilization2(double *, double *, double *, double *, int, int);
template void gpuAddStabilization2(float *, float *, float *, float *, int, int);


template <typename T>   
__global__ void gpuTemplateAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu)
{	
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < M) {
        int i = tid%numPoints;   
        int j = (tid-i)/numPoints;
        for (int k=0; k<ncu; k++) {
            int nm = k * ncu + j;
            int nk = k * numPoints + i;
            fg[tid] += tau[nm] * (ug1[nk] - ug2[nk]);
        }
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAddStabilization3(T *fg, T *ug1, T *ug2, T *tau, int M, int numPoints, int ncu)
{
    int blockDim = 256;
    int gridDim = (M + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateAddStabilization3<<<gridDim, blockDim>>>(fg, ug1, ug2, tau, M, numPoints, ncu);
}

template void gpuAddStabilization3(double *, double *, double *, double *, int, int, int);
template void gpuAddStabilization3(float *, float *, float *, float *, int, int, int);


#endif

