#ifndef __GPUAPP
#define __GPUAPP

// If getting unexpected NaNs when compiling with Enzyme, consider running the script
//      sed -i 's/sqrt/halfpow/g' *.cu
// in the app folder. At the time of this pull request, sqrt(x) will return 0 when compiled with Enzyme, but powf(x, 0.5) works as expected.

template <typename T> __device__ T halfpow(T x){
        return powf(x, 0.5);
}
#ifdef _ENZYME

///// Flux and source
void __device__ __enzyme_fwddiffElement(void*, 
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *, double *,
                int, double *, double *,
                int, double *,
                int, double *,
                int, double,
                int, int, 
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int);

void __device__ __enzyme_fwddiffSource(void*, 
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *,
                int, double,
                int, int, 
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int);

///// ODG type functions: assumes no dependence on dodg
///// TODO: might be okay to just make sure dodg is set to 0 here
void __device__ __enzyme_fwddiffElementField(void*, 
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *,
                int, double,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int);

/// Ubou and Fbou
void __device__ __enzyme_fwddiffgpuUbou(void*, 
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *,
                int, double *,
                int, double *,
                int, double *,
                int, double *,
                int, double,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int);

void __device__ __enzyme_fwddiffgpuFbou(void*, 
                int, double *, double *,
                int, double *,
                int, double *, double *,
                int, double *, double *,
                int, double *, double *,
                int, double *, double *,
                int, double *,
                int, double *,
                int, double *,
                int, double *,
                int, double,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int,
                int, int);     

// //// Fhat and stab
// TODO: Needed for custom Fhat and Stab; 
// void __device__ __enzyme_fwddiff(void*, 
//         int, double *, double *,
//         int, double *,
//         int, double *, double *,
//         int, double *, double *,
//         int, double *,
//         int, double *,
//         int, double *, double *,
//         int, double *, double *,
//         int, double *, double *, 
//         int, double *,
//         int, double *,
//         int, double *,
//         int, double *, 
//         int, double,
//         int, int,
//         int, int,
//         int, int,
//         int, int,
//         int, int,
//         int, int,
//         int, int,
//         int, int);

int __device__ enzyme_const;
int __device__ enzyme_dup;
#endif

#include <stdio.h>
#include <math.h>

#include "gpuFlux.cu"
#include "gpuSource.cu"
#include "gpuSourcew.cu"
#include "gpuEoS.cu"
#include "gpuEoSdu.cu"
#include "gpuEoSdw.cu"
#include "gpuTdfunc.cu"
#include "gpuAvfield.cu"
#include "gpuOutput.cu"
#include "gpuFbou.cu"
#include "gpuUbou.cu"
#include "gpuFhat.cu"
#include "gpuUhat.cu"
#include "gpuStab.cu"
#include "gpuInitu.cu"
#include "gpuInitq.cu"
#include "gpuInitudg.cu"
#include "gpuInitwdg.cu"
#include "gpuInitodg.cu"

#ifdef _ENZYME

///// Flux
template <typename T> __global__ void kernelgpuGradFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
    __enzyme_fwddiffElement((void*)devicegpuFlux<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg, dudg,
                enzyme_dup, odg, dodg,
                enzyme_dup, wdg, dwdg,
                enzyme_const, uinf,
                enzyme_const, param,
                enzyme_const, time,
                enzyme_const, modelnumber, 
                enzyme_const, ng,
                enzyme_const, nc,
                enzyme_const, ncu,
                enzyme_const, nd,
                enzyme_const, ncx,
                enzyme_const, nco,
                enzyme_const, ncw);            
}

template <typename T> void gpuFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuGradFluxEnzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, dodg, wdg, dwdg, uinf, param, 
                                                        time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFluxEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);

///// Source
template <typename T> __global__ void kernelgpuGradSourceEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
    __enzyme_fwddiffSource((void*)devicegpuSource<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg, dudg,
                enzyme_const, odg,
                enzyme_dup, wdg, dwdg,
                enzyme_const, uinf,
                enzyme_const, param,
                enzyme_const, time,
                enzyme_const, modelnumber, 
                enzyme_const, ng,
                enzyme_const, nc,
                enzyme_const, ncu,
                enzyme_const, nd,
                enzyme_const, ncx,
                enzyme_const, nco,
                enzyme_const, ncw);            
}

template <typename T> void gpuSourceEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuGradSourceEnzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uinf, param, 
        time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuSourceEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);
   
// Ubou and Fbou: code is generated because it depends on number of boundaries

template <typename T> __global__ void kernelgpuGradAvfieldEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)
{    
    __enzyme_fwddiffElementField((void*)devicegpuAvfield<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg, dudg,
                enzyme_const, odg,
                enzyme_dup, wdg, dwdg,
                enzyme_const, uinf,
                enzyme_const, param,
                enzyme_const, time,
                enzyme_const, modelnumber,
                enzyme_const, ng,
                enzyme_const, nc,
                enzyme_const, ncu,
                enzyme_const, nd,
                enzyme_const, ncx,
                enzyme_const, nco,
                enzyme_const, ncw,
                enzyme_const, nce,
                enzyme_const, npe,
                enzyme_const, ne);            
}

template <typename T> void gpuAvfieldEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuGradAvfieldEnzyme<<<gridDim, blockDim>>>(f,df, xg, udg, dudg, odg, wdg, dwdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, nce, npe, ne);
}

template void gpuAvfieldEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);

#endif


#endif

