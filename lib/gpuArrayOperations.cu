#ifndef __GPUARRAYOPERATIONS
#define __GPUARRAYOPERATIONS

template <typename T>
__global__ void gpuPrintArray2D(T* a, int m, int n)
{        
    for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++)
            printf("%g   ", a[j*m+i]);                
        printf("\n");
    }
    printf("\n");
}

template <typename T> void gpuPrint2DArray(T* a, int m, int n)
{
    gpuPrintArray2D<<<1, 1>>>(a, m, n);
}

template <typename T>
__global__ void gpuPrintArray3D(T* a, int m, int n, int p)
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

template <typename T> void gpuPrint3DArray(T* a, int m, int n, int p)
{
    gpuPrintArray3D<<<1, 1>>>(a, m, n, p);
}


template <typename T>
__global__ void gpuTemplateGetArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[ind[tid]];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGetArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateGetArrayAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuTemplatePutArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] = x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuPutArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplatePutArrayAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}


template <typename T>
__global__ void gpuTemplateArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += a*x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXPYAtIndex<<<gridDim, blockDim>>>(y, x, a, ind, n);
}

template <typename T>
__global__ void gpuTemplateArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayPlusXAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuTemplateArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] -= x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayMinusXAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuTemplateArraySetValueAtIndex(T *y, T a, int n)
{    
    y[n] = a;          
}

template <typename T> void gpuArraySetValueAtIndex(T *y, T a, int n)
{            
    gpuTemplateArraySetValueAtIndex<<<1, 1>>>(y, a, n);
}

template <typename T>
__global__ void gpuTemplateArraySetValue(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySetValue(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArraySetValue<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuTemplateArrayMultiplyScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMultiplyScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayMultiplyScalar<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuTemplateArrayAddScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] += a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAddScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAddScalar<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuTemplateArrayCopy(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCopy(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayCopy<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayMinus(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = -x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMinus(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayMinus<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayAbs(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = fabs(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAbs(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAbs<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArraySqrt(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sqrt(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySqrt(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArraySqrt<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArraySin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArraySin<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayCos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayCos<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayTan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayTan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayTan<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayAsin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAsin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAsin<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayAcos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAcos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAcos<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayAtan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAtan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAtan<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArraySinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArraySinh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayCosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayCosh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayTanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayTanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayTanh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayAsinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAsinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAsinh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuTemplateArrayAcosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAcosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAcosh<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayAtanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAtanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAtanh<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayExp(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = exp(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayExp(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayExp<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayLog(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = log(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayLog(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayLog<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayCeil(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = ceil(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCeil(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayCeil<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayFloor(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = floor(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayFloor(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayFloor<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayErf(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erf(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayErf(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayErf<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayErfc(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erfc(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayErfc(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayErfc<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArraySquare(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid]*x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySquare(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArraySquare<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuTemplateArrayPower(T *y, T *x, int p, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        for (int j=1; j<p; j++)
            y[tid] = y[tid]*x[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayPower(T *y, T *x, int p, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayPower<<<gridDim, blockDim>>>(y, x, p, n);
}


template <typename T>
__global__ void gpuTemplateArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] = a*C[tid+n*tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void gpuArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayMultiplyScalarDiagonal<<<gridDim, blockDim>>>(C, a, n);
}


template <typename T>
__global__ void gpuTemplateArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] += a*x[tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void gpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAddVectorToDiagonal<<<gridDim, blockDim>>>(C, x, a, n);
}


template <typename T>
__global__ void gpuTemplateArrayRowAverage(T *y, T *x, int m, int n)
{    
    int j = threadIdx.x + blockIdx.x * blockDim.x;    
    while (j < n) {        
        T avg = 0;
        int i;
        for (i=0; i<m; i++)
            avg = avg + x[i + m*j];
        avg = avg/((T) m);
        for (i=0; i<m; i++)
            y[i + m*j] = avg;         
        j += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuArrayRowAverage(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayRowAverage<<<gridDim, blockDim>>>(y, x, m, n);
}


template <typename T>
__global__ void gpuTemplateArrayAXPB(T *y, T *x, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*x[tid]+b;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPB(T *y, T *x, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXPB<<<gridDim, blockDim>>>(y, x, a, b, n);
}


template <typename T>
__global__ void gpuTemplateArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        z[tid] = a*x[tid]+b*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXPBY<<<gridDim, blockDim>>>(z, x, y, a, b, n);
}

template <typename T>
__global__ void gpuTemplateArrayAXY(T *s, T *x, T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXY(T *s, T *x, T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXY<<<gridDim, blockDim>>>(s, x, y, a, n);
}

template <typename T>
__global__ void gpuTemplateArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid]*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXYZ<<<gridDim, blockDim>>>(s, x, y, z, a, n);
}

template <typename T>
__global__ void gpuTemplateArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid] + b*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAXYPBZ<<<gridDim, blockDim>>>(s, x, y, z, a, b, n);
}

template <typename T>
__global__ void gpuTemplateArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid] + b*y[tid] + c*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAdd3Vectors<<<gridDim, blockDim>>>(s, x, y, z, a, b, c, n);
}

template <typename T>
__global__ void  gpuTemplateArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        a[tid] = b[tid] + c[tid] + d[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayAdd3Vector<<<gridDim, blockDim>>>(a, b, c, d, n);
}

template <typename T>
__global__ void gpuTemplateArrayExtract(T *un, T *u, int I, int J, int M, int N,  
        int i1, int j1, int k1, int ni)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+I*J*k];      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayExtract<<<gridDim, blockDim>>>(un, u, I, J, M, N, i1, j1, k1, ni);
}

template <typename T>
__global__ void gpuTemplateArrayInsert(T *u, T *un, int I, int J, int M, int N,  
        int i1, int j1, int k1, int ni)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+I*J*k] = un[idx];      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayInsert<<<gridDim, blockDim>>>(u, un, I, J, M, N, i1, j1, k1, ni);
}

template <typename T> 
__global__ void gpuTemplateArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K, int N, int Q)
{        
    // static shared memory
    __shared__ T Ashared[256];

    if (threadIdx.x<Q)
    {
      // load data from global memory to shared memory
      Ashared[threadIdx.x] = A[threadIdx.x];
    }

    // thread synchronization
    __syncthreads();

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {        
        int i = idx%I;
        int j = (idx-i)/I;                
        int m = K*j;
        C[idx] = 0.0;
        for (int k=0; k<K; k++)
            C[idx] += Ashared[i+I*k]*B[k+m];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K)
{        
    // C[I*J] = A[I*K] x B[K*J]
    int N = I*J;    
    int Q = I*K;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayGemmSharedMem<<<gridDim, blockDim>>>(C, A, B, I, J, K, N, Q);
}

template <typename T> 
__global__ void gpuTemplateArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        int a = i+Q*s;
        int b = K*j+P*s;
        C[idx] = 0.0;
        for (int k=0; k<K; k++)
            C[idx] += A[a+I*k]*B[k+b];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayGemmBatch<<<gridDim, blockDim>>>(C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void gpuTemplateArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        int a = i+Q*s;
        int b = K*j+P*s;
        for (int k=0; k<K; k++)
            C[idx] += A[a+I*k]*B[k+b];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayGemmBatch1<<<gridDim, blockDim>>>(C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void gpuTemplateArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nent) {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++)
            ucg[i] += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = ucg[i]/((T) nelem);        
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayDG2CG<<<gridDim, blockDim>>>(ucg, udg, cgent2dgent, rowent2elem, nent);
}

template <typename T> 
__global__ void gpuTemplateArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nent) {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            for (int j=0; j<npe; j++)
                ucg[i] += udg[j+npe*e]; 
        }
        ucg[i] = ucg[i]/((T) (nelem*npe));
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayDG2CG2<<<gridDim, blockDim>>>(ucg, udg, colent2elem, rowent2elem, nent, npe);
}

template <typename T> __global__ void gpuTemplateArrayInverseMatrix11(T *A, int N)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        A[i] = 1.0/A[i];    
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayInverseMatrix11(T *A, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayInverseMatrix11<<<gridDim, blockDim>>>(A, N);
}

template <typename T> __global__ void gpuTemplateArrayInverseMatrix22(T *A, int N)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        double a11 = A[i + N*0];
        double a21 = A[i + N*1];
        double a12 = A[i + N*2];
        double a22 = A[i + N*3];
        double detA = (a11*a22- a12*a21);
      
        A[i + N*0] = a22/detA;
        A[i + N*1] = -a21/detA;
        A[i + N*2] = -a12/detA;
        A[i + N*3] = a11/detA;
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayInverseMatrix22(T *A, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayInverseMatrix22<<<gridDim, blockDim>>>(A, N);
}

template <typename T> __global__ void gpuTemplateArrayInverseMatrix33(T *A, int N)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        double a11 = A[i + N*0];
        double a21 = A[i + N*1];
        double a31 = A[i + N*2];
        double a12 = A[i + N*3];
        double a22 = A[i + N*4];
        double a32 = A[i + N*5];
        double a13 = A[i + N*6];
        double a23 = A[i + N*7];
        double a33 = A[i + N*8];        
        double detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
      
        A[i + N*0] = (a22*a33 - a23*a32)/detA;
        A[i + N*1] = (a23*a31 - a21*a33)/detA;
        A[i + N*2] = (a21*a32 - a22*a31)/detA;
        A[i + N*3] = (a13*a32 - a12*a33)/detA;
        A[i + N*4] = (a11*a33 - a13*a31)/detA;
        A[i + N*5] = (a12*a31 - a11*a32)/detA;
        A[i + N*6] = (a12*a23 - a13*a22)/detA;
        A[i + N*7] = (a13*a21 - a11*a23)/detA;
        A[i + N*8] = (a11*a22 - a12*a21)/detA;        
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayInverseMatrix33(T *A, int N)
{
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayInverseMatrix33<<<gridDim, blockDim>>>(A, N);
}

template <typename T> __global__ void gpuTemplateArrayMatrixMultiplication(T *C, T *A, T *B, 
        int S, int I, int J, int K, int M, int N)
{        
    // C[S*I*J] = A[S*I*K] x B[S*K*J]
    //int M = I*J;
    //int N = M*S;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M; //   [1, I*J]        
        int i = l%I;   //   [1, I]
        int j = (l-i)/I; // [1, J]
        int s = (idx-l)/M;//[1, S]  
        C[s + S*i + S*I*j] = 0.0;
        for (int k=0; k<K; k++)
            C[s + S*i + S*I*j] += A[s + S*i + S*I*k]*B[s + S*k + S*K*j];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayMatrixMultiplication(T *C, T *A, T *B, int S, int I, int J, int K)
{        
    int M = I*J;
    int N = M*S;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayMatrixMultiplication<<<gridDim, blockDim>>>(C, A, B, S, I, J, K, M, N);
}

template <typename T> __global__ void gpuTemplateArrayEosInverseMatrix11(T *A, int N)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        A[i] = 1.0/A[i];    
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayEosInverseMatrix11(T *A, int npe, int ncw, int ne)
{
    int N = npe*ne;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayEosInverseMatrix11<<<gridDim, blockDim>>>(A, int N);
}

template <typename T> __global__ void gpuTemplateArrayEosInverseMatrix22(T *A, int N, int M, int npe)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        int j = i%npe;    // [1, npe]
        int k = (i-j)/npe; //[1, ne]
        double a11 = A[j + npe*0 + M*k];
        double a21 = A[j + npe*1 + M*k];
        double a12 = A[j + npe*2 + M*k];
        double a22 = A[j + npe*3 + M*k];
        double detA = (a11*a22- a12*a21);      
        A[j + npe*0 + M*k] = a22/detA;
        A[j + npe*1 + M*k] = -a21/detA;
        A[j + npe*2 + M*k] = -a12/detA;
        A[j + npe*3 + M*k] = a11/detA;
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayEosInverseMatrix22(T *A, int npe, int ncw, int ne)
{
    int N = npe*ne;
    int M = npe*ncw*ncw;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayEosInverseMatrix22<<<gridDim, blockDim>>>(A, N, M, npe);
}

template <typename T> __global__ void gpuTemplateArrayEosInverseMatrix33(T *A, int N, int M, int npe)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < N) {
        int j = i%npe;    // [1, npe]
        int k = (i-j)/npe; //[1, ne]
        double a11 = A[j + npe*0 + M*k];
        double a21 = A[j + npe*1 + M*k];
        double a31 = A[j + npe*2 + M*k];
        double a12 = A[j + npe*3 + M*k];
        double a22 = A[j + npe*4 + M*k];
        double a32 = A[j + npe*5 + M*k];
        double a13 = A[j + npe*6 + M*k];
        double a23 = A[j + npe*7 + M*k];
        double a33 = A[j + npe*8 + M*k];        
        double detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);
      
        A[j + npe*0 + M*k] = (a22*a33 - a23*a32)/detA;
        A[j + npe*1 + M*k] = (a23*a31 - a21*a33)/detA;
        A[j + npe*2 + M*k] = (a21*a32 - a22*a31)/detA;
        A[j + npe*3 + M*k] = (a13*a32 - a12*a33)/detA;
        A[j + npe*4 + M*k] = (a11*a33 - a13*a31)/detA;
        A[j + npe*5 + M*k] = (a12*a31 - a11*a32)/detA;
        A[j + npe*6 + M*k] = (a12*a23 - a13*a22)/detA;
        A[j + npe*7 + M*k] = (a13*a21 - a11*a23)/detA;
        A[j + npe*8 + M*k] = (a11*a22 - a12*a21)/detA;            
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayEosInverseMatrix33(T *A, int npe, int ncw, int ne)
{
    int N = npe*ne;
    int M = npe*ncw*ncw;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayEosInverseMatrix33<<<gridDim, blockDim>>>(A, N, M, npe);
}

template <typename T> __global__ void gpuTemplateArrayEosMatrixMultiplication(T *C, T *A, T *B, 
        int npe, int ncw, int ncu, int N, int K, int P, int Q)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int j = idx%npe;    // [1, npe]
        int k = (i-j)/npe;  // [1, ne]        
        for (int b=0; b<ncu; b++)
          for (int a=0; a<ncw; a++) {
            C[j + npe*a + K*b + P*k] = 0.0;
            for (int m=0; m<ncw; m++)
              C[j + npe*a + K*b + P*k] += A[j + npe*a + K*m + Q*k]*B[j + npe*m + K*b + P*k];
          }        
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void gpuArrayEosMatrixMultiplication(T *C, T *A, T *B, int npe, int ncw, int ne, int ncu)
{        
    int N = npe*ne;
    int K = npe*ncw;
    int P = K*ncu;
    int Q = K*ncw;    
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateArrayEosMatrixMultiplication<<<gridDim, blockDim>>>(C, A, B, npe, ncw, ncu, N, K, P, Q);
}

template void gpuPrint2DArray(double*, int, int);
template void gpuPrint3DArray(double*, int, int, int);
template void gpuArraySetValue(double*, double, int);
template void gpuArraySetValueAtIndex(double*, double, int);
template void gpuArrayAddScalar(double*, double, int);
template void gpuArrayMultiplyScalar(double*, double, int);

template void gpuGetArrayAtIndex(double*, double*, int*, int);
template void gpuPutArrayAtIndex(double*, double*, int*, int);
template void gpuArrayAXPYAtIndex(double*, double*, double, int*, int);
template void gpuArrayPlusXAtIndex(double*, double*, int*, int);
template void gpuArrayMinusXAtIndex(double*, double*, int*, int);

template void gpuArrayCopy(double*, double*, int);
template void gpuArrayMinus(double*, double*, int);
template void gpuArrayAbs(double*, double*, int);
template void gpuArraySqrt(double*, double*, int);
template void gpuArraySin(double*, double*, int);
template void gpuArrayCos(double*, double*, int);
template void gpuArrayTan(double*, double*, int);
template void gpuArrayAsin(double*, double*, int);
template void gpuArrayAcos(double*, double*, int);
template void gpuArrayAtan(double*, double*, int);
template void gpuArraySinh(double*, double*, int);
template void gpuArrayCosh(double*, double*, int);
template void gpuArrayTanh(double*, double*, int);
template void gpuArrayAsinh(double*, double*, int);
template void gpuArrayAcosh(double*, double*, int);
template void gpuArrayAtanh(double*, double*, int);
template void gpuArrayExp(double*, double*, int);
template void gpuArrayLog(double*, double*, int);
template void gpuArrayCeil(double*, double*, int);
template void gpuArrayFloor(double*, double*, int);
template void gpuArrayErf(double*, double*, int);
template void gpuArrayErfc(double*, double*, int);
template void gpuArraySquare(double*, double*, int);
template void gpuArrayPower(double*, double*, int, int);

template void gpuArrayMultiplyScalarDiagonal(double*, double, int);
template void gpuArrayAddVectorToDiagonal(double*, double*, double, int);
template void gpuArrayRowAverage(double*, double*, int, int);
template void gpuArrayAXPB(double*, double*, double, double, int);
template void gpuArrayAXPBY(double*, double*, double*, double, double, int);
template void gpuArrayAXY(double*, double*, double*, double, int);
template void gpuArrayAXYZ(double*, double*, double*, double*, double, int);
template void gpuArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void gpuArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void gpuArrayAdd3Vector(double*, double*, double*, double*, int);
template void gpuArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void gpuArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void gpuArrayGemmSharedMem(double*, double*, double*, int, int, int);
template void gpuArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void gpuArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void gpuArrayDG2CG(double*, double*, int*, int*, int);
template void gpuArrayDG2CG2(double*, double*, int*, int*, int, int);

template void gpuPrint2DArray(float*, int, int);
template void gpuPrint3DArray(float*, int, int, int);
template void gpuArraySetValue(float*, float, int);
template void gpuArraySetValueAtIndex(float*, float, int);
template void gpuArrayAddScalar(float*, float, int);
template void gpuArrayMultiplyScalar(float*, float, int);

template void gpuGetArrayAtIndex(float*, float*, int*, int);
template void gpuPutArrayAtIndex(float*, float*, int*, int);
template void gpuArrayAXPYAtIndex(float*, float*, float, int*, int);
template void gpuArrayPlusXAtIndex(float*, float*, int*, int);
template void gpuArrayMinusXAtIndex(float*, float*, int*, int);

template void gpuArrayCopy(float*, float*, int);
template void gpuArrayMinus(float*, float*, int);
template void gpuArrayAbs(float*, float*, int);
template void gpuArraySqrt(float*, float*, int);
template void gpuArraySin(float*, float*, int);
template void gpuArrayCos(float*, float*, int);
template void gpuArrayTan(float*, float*, int);
template void gpuArrayAsin(float*, float*, int);
template void gpuArrayAcos(float*, float*, int);
template void gpuArrayAtan(float*, float*, int);
template void gpuArraySinh(float*, float*, int);
template void gpuArrayCosh(float*, float*, int);
template void gpuArrayTanh(float*, float*, int);
template void gpuArrayAsinh(float*, float*, int);
template void gpuArrayAcosh(float*, float*, int);
template void gpuArrayAtanh(float*, float*, int);
template void gpuArrayExp(float*, float*, int);
template void gpuArrayLog(float*, float*, int);
template void gpuArrayCeil(float*, float*, int);
template void gpuArrayFloor(float*, float*, int);
template void gpuArrayErf(float*, float*, int);
template void gpuArrayErfc(float*, float*, int);
template void gpuArraySquare(float*, float*, int);
template void gpuArrayPower(float*, float*, int, int);

template void gpuArrayMultiplyScalarDiagonal(float*, float, int);
template void gpuArrayAddVectorToDiagonal(float*, float*, float, int);
template void gpuArrayRowAverage(float*, float*, int, int);
template void gpuArrayAXPB(float*, float*, float, float, int);
template void gpuArrayAXPBY(float*, float*, float*, float, float, int);
template void gpuArrayAXY(float*, float*, float*, float, int);
template void gpuArrayAXYZ(float*, float*, float*, float*, float, int);
template void gpuArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void gpuArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void gpuArrayAdd3Vector(float*, float*, float*, float*, int);
template void gpuArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void gpuArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void gpuArrayGemmSharedMem(float*, float*, float*, int, int, int);
template void gpuArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void gpuArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void gpuArrayDG2CG(float*, float*, int*, int*, int);
template void gpuArrayDG2CG2(float*, float*, int*, int*, int, int);

template void gpuArrayInverseMatrix11(double*, int);
template void gpuArrayInverseMatrix11(float*, int);
template void gpuArrayInverseMatrix22(double*, int);
template void gpuArrayInverseMatrix22(float*, int);
template void gpuArrayInverseMatrix33(double*, int);
template void gpuArrayInverseMatrix33(float*, int);
template void gpuArrayMatrixMultiplication(double*, double*, double*, int, int, int, int);
template void gpuArrayMatrixMultiplication(float*, float*, float*, int, int, int, int);

template void gpuArrayEosInverseMatrix11(double*, int, int, int);
template void gpuArrayEosInverseMatrix11(float*, int, int, int);
template void gpuArrayEosInverseMatrix22(double*, int, int, int);
template void gpuArrayEosInverseMatrix22(float*, int, int, int);
template void gpuArrayEosInverseMatrix33(double*, int, int, int);
template void gpuArrayEosInverseMatrix33(float*, int, int, int);
template void gpuArrayEosMatrixMultiplication(double*, double*, double*, int, int, int, int);
template void gpuArrayEosMatrixMultiplication(float*, float*, float*, int, int, int, int);

#endif


