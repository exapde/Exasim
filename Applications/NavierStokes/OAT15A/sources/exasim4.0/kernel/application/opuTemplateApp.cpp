#ifndef __OPUAPP
#define __OPUAPP

#include <math.h>

#include "opuFluxApplicationName.cpp"
#include "opuSourceApplicationName.cpp"
#include "opuTdfuncApplicationName.cpp"
#include "opuFbouApplicationName.cpp"
#include "opuUbouApplicationName.cpp"
#include "opuUboutdepApplicationName.cpp"
#include "opuAVfieldApplicationName.cpp"
#include "opuStabApplicationName.cpp"

template <typename T> void opuFlux(T *f, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{            
    opuFluxApplicationName(f, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);        
}
template void opuFlux(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void opuFlux(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void opuSource(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuSourceApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void opuSource(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void opuSource(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T> void opuTdfunc(T *s, T *xg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuTdfuncApplicationName(s, xg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco);   
}
template void opuTdfunc(double *, double *, double *, double *, double *, double, 
            int, int, int, int, int, int, int);
template void opuTdfunc(float *, float *, float *, float *, float *, float, 
            int, int, int, int, int, int, int);

template <typename T>   
void opuFbou(T *fh, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuFbouApplicationName(fh, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void opuFbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuFbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T>   
void opuUbou(T *ub, T *xg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuUbouApplicationName(ub, xg, udg, uhg, odg, nl, tau, uinf, param, time, 
                 ib, ng, nc, ncu, nd, ncx, nco);
}
template void opuUbou(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuUbou(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T> 
void opuAVfield(T *f, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne)
{
	opuAVfieldApplicationName(f, xdg, udg, odg, param, time, ng, nc, ncu, nd, ncx, nco, npe, ne);
}
template void opuAVfield(double *, double *, double *, double *, double *, double, 
        int, int, int, int, int, int, int, int);
template void opuAVfield(float *, float *, float *, float *, float *, float, 
        int, int, int, int, int, int, int, int);

template <typename T>   
void opuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nlg, T *param, 
        T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuStabApplicationName(fh, xg, uhg, udg1, udg2, odg1, odg2, nlg, param, time, 
                 ng, nc, ncu, nd, ncx, nco);
}
template void opuStab(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double, int, int, int, int, int, int, int);
template void opuStab(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float, int, int, int, int, int, int, int);

template <typename T>   
void opuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt,
        int ib, int jb, int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname)
{        
    opuUboutdepApplicationName(ub, xg, udg, uhg, odg, udg0, udg1, udg2, uhg0, uhg1, uhg2, 
            nl, tau, uinf, param, time, dt, ib, jb, tstage, ng, nc, ncu, nd, ncx, nco);
}
template void opuUboutdep(double *, double *, double *, double *, double *, double *,
        double *, double *, double *, double *, double *, double *, double *, double *, 
        double *, double, double, int, int, int, int, int, int, int, int, int, int);
template void opuUboutdep(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float *, float *, float *, float *, float *, 
        float *, float, float, int, int, int, int, int, int, int, int, int, int);

#endif

