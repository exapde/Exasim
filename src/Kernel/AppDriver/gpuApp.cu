#ifndef __GPUAPP
#define __GPUAPP

#ifdef _ENZYME
// template <typename... Args>
// void __device__ __enzyme_autodiff(void*, Args... args);

///// Flux and source
void __device__ __enzyme_fwddiffElement(void*, 
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
                int, double *,
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
//#include "gpuEoS.cu"
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
// template <typename... Args>
// // void __device__ __enzyme_autodiff(void*, Args... args);

// ///// Flux and source
// void __device__ __enzyme_fwddiffElement(void*, 
//                 int, double *, double *,
//                 int, double *,
//                 int, double *, double *,
//                 int, double *,
//                 int, double *, double *,
//                 int, double *,
//                 int, double *,
//                 int, double,
//                 int, int, 
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int);

// /// Ubou and Fbou
// void __device__ __enzyme_fwddiffFace(void*, 
//                 int, double *, double *,
//                 int, double *,
//                 int, double *, double *,
//                 int, double *,
//                 int, double *, double *,
//                 int, double *,
//                 int, double *,
//                 int, double *,
//                 int, double *,
//                 int, double *,
//                 int, double,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int,
//                 int, int);   

// // //// Fhat and stab
// // void __device__ __enzyme_fwddiff(void*, 
// //         int, double *, double *,
// //         int, double *,
// //         int, double *, double *,
// //         int, double *, double *,
// //         int, double *,
// //         int, double *,
// //         int, double *, double *,
// //         int, double *, double *,
// //         int, double *, double *, 
// //         int, double *,
// //         int, double *,
// //         int, double *,
// //         int, double *, 
// //         int, double,
// //         int, int,
// //         int, int,
// //         int, int,
// //         int, int,
// //         int, int,
// //         int, int,
// //         int, int,
// //         int, int);

// int __device__ enzyme_const;
// int __device__ enzyme_dup;

///// Flux
template <typename T> __global__ void kernelgpuGradFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
//     printf("Calling gpuFLUXEnzyme\n");
    __enzyme_fwddiffElement((void*)devicegpuFlux<T>, 
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

template <typename T> void gpuFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
        // printf("do i get inside gpufluxenzyme\n");
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
        // printf("before kernel\n");
	kernelgpuGradFluxEnzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uinf, param, 
        time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
        // printf("after kernel\n");
}

template void gpuFluxEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);

///// Source
template <typename T> __global__ void kernelgpuGradSourceEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
//     printf("Calling gpuFLUXEnzyme\n");
    __enzyme_fwddiffElement((void*)devicegpuSource<T>, 
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
        // printf("after kernel\n");
}

template void gpuSourceEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);

// ///// Ubou
// template <typename T> __global__ void kernelgpuGradUbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, 
//         T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
// {    
//         if (ib == 1){
//                 __enzyme_fwddiffFace((void*)devicegpuUbou1<T>, 
//                                 enzyme_dup, f, df,
//                                 enzyme_const, xg,
//                                 enzyme_dup, udg, dudg,
//                                 enzyme_const, odg,
//                                 enzyme_dup, wdg, dwdg,
//                                 enzyme_const, uhg,
//                                 enzyme_const, nlg,
//                                 enzyme_const, tau,
//                                 enzyme_const, uinf,
//                                 enzyme_const, param,
//                                 enzyme_const, time,
//                                 enzyme_const, modelnumber, 
//                                 enzyme_const, ng,
//                                 enzyme_const, nc,
//                                 enzyme_const, ncu,
//                                 enzyme_const, nd,
//                                 enzyme_const, ncx,
//                                 enzyme_const, nco,
//                                 enzyme_const, ncw);    
//         }
        
// }

// template <typename T> void gpuUbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, 
//         T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
// {
// 	int blockDim = 256;
// 	int gridDim = (ng + blockDim - 1) / blockDim;
// 	gridDim = (gridDim>1024)? 1024 : gridDim;
// 	if (ib == 1)
// 		kernelgpuGradUbouEnzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
// }

// template void gpuUbouEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
//         double, int, int, int, int, int, int, int, int, int);

///// Fbou
// template <typename T> __global__ void kernelgpuGradFbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, 
//         T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
// {    
//         if (ib == 1){
//                 __enzyme_fwddiffFace((void*)devicegpuFbou1<T>, 
//                                 enzyme_dup, f, df,
//                                 enzyme_const, xg,
//                                 enzyme_dup, udg, dudg,
//                                 enzyme_const, odg,
//                                 enzyme_dup, wdg, dwdg,
//                                 enzyme_const, uhg,
//                                 enzyme_const, nlg,
//                                 enzyme_const, tau,
//                                 enzyme_const, uinf,
//                                 enzyme_const, param,
//                                 enzyme_const, time,
//                                 enzyme_const, modelnumber, 
//                                 enzyme_const, ng,
//                                 enzyme_const, nc,
//                                 enzyme_const, ncu,
//                                 enzyme_const, nd,
//                                 enzyme_const, ncx,
//                                 enzyme_const, nco,
//                                 enzyme_const, ncw);    
//         }
        
// }

// template <typename T> void gpuFbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
//         T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
// {
// 	int blockDim = 256;
// 	int gridDim = (ng + blockDim - 1) / blockDim;
// 	gridDim = (gridDim>1024)? 1024 : gridDim;
// 	if (ib == 1)
// 		kernelgpuGradFbouEnzyme<<<gridDim, blockDim>>>(f, df, xg, udg, dudg, odg, wdg, dwdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
// }

// template void gpuFbouEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double  *, 
//         double, int, int, int, int, int, int, int, int, int);
#endif


#endif

