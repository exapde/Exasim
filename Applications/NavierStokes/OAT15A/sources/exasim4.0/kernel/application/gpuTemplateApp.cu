#ifndef __GPUAPP
#define __GPUAPP

#include "gpuFluxApplicationName.cu"
#include "gpuSourceApplicationName.cu"
#include "gpuTdfuncApplicationName.cu"
#include "gpuFbouApplicationName.cu"
#include "gpuUbouApplicationName.cu"
#include "gpuUboutdepApplicationName.cu"
#include "gpuAVfieldApplicationName.cu"
#include "gpuStabApplicationName.cu"

template <typename T> void gpuFlux(T *f, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{            
    gpuFluxApplicationName(f, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);        
}
template void gpuFlux(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void gpuFlux(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void gpuSource(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuSourceApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void gpuSource(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void gpuSource(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void gpuTdfunc(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuTdfuncApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void gpuTdfunc(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void gpuTdfunc(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T>   
void gpuFbou(T *fh, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuFbouApplicationName(fh, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void gpuFbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuFbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T>   
void gpuUbou(T *ub, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuUbouApplicationName(ub, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void gpuUbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuUbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T> 
void gpuAVfield(T *f, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne)
{
	gpuAVfieldApplicationName(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);
}
template void gpuAVfield(double *, double *, double *, double *, double *, double, 
        int, int, int, int, int, int, int, int);
template void gpuAVfield(float *, float *, float *, float *, float *, float, 
        int, int, int, int, int, int, int, int);

template <typename T>   
void gpuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nlg, T *param, 
        T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuStabApplicationName(fh, xg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, 
                 ng, nc, ncu, nd, ncx, nco);
}
template void gpuStab(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int);
template void gpuStab(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int);

template <typename T>   
void gpuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt,
        int ib, int jb, int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    gpuUboutdepApplicationName(ub, xg, udg, uhg, odg, udg0, udg1, udg2, uhg0, uhg1, uhg2, 
            nl, tau, uinf, param, time, dt, ib, jb, tstage, ng, nc, ncu, nd, ncx, nco);
}
template void gpuUboutdep(double *, double *, double *, double *, double *, double *,
        double *, double *, double *, double *, double *, double *, double *, double *, 
        double *, double, double, int, int, int, int, int, int, int, int, int, int);
template void gpuUboutdep(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float *, float *, float *, float *, float *, 
        float *, float, float, int, int, int, int, int, int, int, int, int, int);

#endif

