#ifndef __CPUAPP
#define __CPUAPP

#include <math.h>
#include <omp.h>

#include "cpuFluxApplicationName.cpp"
#include "cpuSourceApplicationName.cpp"
#include "cpuTdfuncApplicationName.cpp"
#include "cpuFbouApplicationName.cpp"
#include "cpuUbouApplicationName.cpp"
#include "cpuUboutdepApplicationName.cpp"
#include "cpuAVfieldApplicationName.cpp"
#include "cpuStabApplicationName.cpp"

template <typename T> void cpuFlux(T *f, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{            
    cpuFluxApplicationName(f, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);        
}
template void cpuFlux(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void cpuFlux(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void cpuSource(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuSourceApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void cpuSource(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void cpuSource(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void cpuTdfunc(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuTdfuncApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void cpuTdfunc(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void cpuTdfunc(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T>   
void cpuFbou(T *fh, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuFbouApplicationName(fh, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void cpuFbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuFbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T>   
void cpuUbou(T *ub, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuUbouApplicationName(ub, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void cpuUbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuUbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T> 
void cpuAVfield(T *f, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne)
{
	cpuAVfieldApplicationName(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);
}
template void cpuAVfield(double *, double *, double *, double *, double *, double, 
        int, int, int, int, int, int, int, int);
template void cpuAVfield(float *, float *, float *, float *, float *, float, 
        int, int, int, int, int, int, int, int);

template <typename T>   
void cpuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nlg, T *param, 
        T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuStabApplicationName(fh, xg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, 
                 ng, nc, ncu, nd, ncx, nco);
}
template void cpuStab(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int);
template void cpuStab(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int);

template <typename T>   
void cpuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt,
        int ib, int jb, int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    cpuUboutdepApplicationName(ub, xg, udg, uhg, odg, udg0, udg1, udg2, uhg0, uhg1, uhg2, 
            nl, tau, uinf, param, time, dt, ib, jb, tstage, ng, nc, ncu, nd, ncx, nco);
}
template void cpuUboutdep(double *, double *, double *, double *, double *, double *,
        double *, double *, double *, double *, double *, double *, double *, double *, 
        double *, double, double, int, int, int, int, int, int, int, int, int, int);
template void cpuUboutdep(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float *, float *, float *, float *, float *, 
        float *, float, float, int, int, int, int, int, int, int, int, int, int);


#endif

