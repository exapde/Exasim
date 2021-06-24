#ifndef __GPUELEMFACEGEOM
#define __GPUELEMFACEGEOM

template <typename T>
__global__ void  gpuTemplateElemGeom1D(T *jac, T *Xx, T *Jg, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg[tid];
        Xx[tid] = 1.0;
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuElemGeom1D(T *jac, T *Xx, T *Jg, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateElemGeom1D<<<gridDim, blockDim>>>(jac, Xx, Jg, nga);
}
template void gpuElemGeom1D(double*, double*, double*, int);
template void gpuElemGeom1D(float*, float*, float*, int);


template <typename T>
__global__ void  gpuTemplateElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22,
                 T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg11[tid]*Jg22[tid] - Jg12[tid]*Jg21[tid];
        Xx11[tid] = Jg22[tid];
        Xx21[tid] = -Jg21[tid];
        Xx12[tid] = -Jg12[tid];
        Xx22[tid] = Jg11[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuElemGeom2D(T *jac, T *Xx11, T *Xx12, T *Xx21, T *Xx22,
                 T *Jg11, T *Jg12, T *Jg21, T *Jg22, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateElemGeom2D<<<gridDim, blockDim>>>(jac, Xx11, Xx12, Xx21, Xx22, Jg11, Jg12, Jg21, Jg22, nga);
}
template void gpuElemGeom2D(double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void gpuElemGeom2D(float*, float*, float*, float*, float*, float*, float*, float*, float*, int);


template <typename T>
__global__ void  gpuTemplateElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, 
                T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33,
                T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, 
                T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < nga) {
        jac[tid] = Jg11[tid]*Jg22[tid]*Jg33[tid] - Jg11[tid]*Jg32[tid]*Jg23[tid] +
                 Jg21[tid]*Jg32[tid]*Jg13[tid] - Jg21[tid]*Jg12[tid]*Jg33[tid] +
                 Jg31[tid]*Jg12[tid]*Jg23[tid] - Jg31[tid]*Jg22[tid]*Jg13[tid];
        Xx11[tid] = Jg22[tid]*Jg33[tid] - Jg23[tid]*Jg32[tid];
        Xx21[tid] = Jg23[tid]*Jg31[tid] - Jg21[tid]*Jg33[tid];
        Xx31[tid] = Jg21[tid]*Jg32[tid] - Jg22[tid]*Jg31[tid];
        Xx12[tid] = Jg13[tid]*Jg32[tid] - Jg12[tid]*Jg33[tid];
        Xx22[tid] = Jg11[tid]*Jg33[tid] - Jg13[tid]*Jg31[tid];
        Xx32[tid] = Jg12[tid]*Jg31[tid] - Jg11[tid]*Jg32[tid];
        Xx13[tid] = Jg12[tid]*Jg23[tid] - Jg13[tid]*Jg22[tid];
        Xx23[tid] = Jg13[tid]*Jg21[tid] - Jg11[tid]*Jg23[tid];
        Xx33[tid] = Jg11[tid]*Jg22[tid] - Jg12[tid]*Jg21[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuElemGeom3D(T *jac, T *Xx11, T *Xx12, T *Xx13, T *Xx21, 
                T *Xx22, T *Xx23, T *Xx31, T *Xx32, T *Xx33,
                T *Jg11, T *Jg12, T *Jg13, T *Jg21, T *Jg22, 
                T *Jg23, T *Jg31, T *Jg32, T *Jg33, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateElemGeom3D<<<gridDim, blockDim>>>(jac, Xx11, Xx12, Xx13, Xx21, Xx22, Xx23, Xx31, Xx32, Xx33, 
            Jg11, Jg12, Jg13, Jg21, Jg22, Jg23, Jg31, Jg32, Jg33, nga);
}
template void gpuElemGeom3D(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*,
                            double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void gpuElemGeom3D(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*,
                            float*, float*, float*, float*, float*, float*, float*, float*, float*, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {        
        jacg[i] = 1.0;
        nlg[i] = -1.0;  
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuFaceGeom1D(T *jacg, T *nlg, T *Jg, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateFaceGeom1D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuFaceGeom1D(double*, double*, double*, int);
template void gpuFaceGeom1D(float*, float*, float*, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {
        int j = i+nga;
        jacg[i] = sqrt(Jg[i]*Jg[i] + Jg[j]*Jg[j]);
        nlg[i] = Jg[j]/jacg[i];
        nlg[j] = -Jg[i]/jacg[i];
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuFaceGeom2D(T *jacg, T *nlg, T *Jg, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateFaceGeom2D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuFaceGeom2D(double*, double*, double*, int);
template void gpuFaceGeom2D(float*, float*, float*, int);


template <typename T>
__global__ void  gpuTemplateFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nga) {
        int j = i+nga;
        int k = i+2*nga;
        int m = i+3*nga;
        int n = i+4*nga;
        int p = i+5*nga;
        nlg[i] = Jg[j]*Jg[p] - Jg[k]*Jg[n];
        nlg[j] = Jg[k]*Jg[m] - Jg[i]*Jg[p];
        nlg[k] = Jg[i]*Jg[n] - Jg[j]*Jg[m];
        jacg[i] = sqrt(nlg[i]*nlg[i] + nlg[j]*nlg[j] + nlg[k]*nlg[k]);
        nlg[i] = nlg[i]/jacg[i];
        nlg[j] = nlg[j]/jacg[i];
        nlg[k] = nlg[k]/jacg[i];
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
void gpuFaceGeom3D(T *jacg, T *nlg, T *Jg, int nga)
{
    int blockDim = 256;
    int gridDim = (nga + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuTemplateFaceGeom3D<<<gridDim, blockDim>>>(jacg, nlg, Jg, nga);
}
template void gpuFaceGeom3D(double*, double*, double*, int);
template void gpuFaceGeom3D(float*, float*, float*, int);

#endif

