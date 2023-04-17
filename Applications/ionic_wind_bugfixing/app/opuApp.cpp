#ifndef __OPUAPP
#define __OPUAPP
#include <stdio.h>
#include <math.h>

#ifdef _MUTATIONPP
#include <mutation++.h>
#include "mutationUtils.cpp"
#endif

#include "opuFlux.cpp"
#include "opuSource.cpp"
#include "opuSourcew.cpp"
#include "opuEoS.cpp"
#include "opuEoSdu.cpp"
#include "opuEoSdw.cpp"
#include "opuTdfunc.cpp"
#include "opuAvfield.cpp"
#include "opuOutput.cpp"
#include "opuFbou.cpp"
#include "opuUbou.cpp"
#include "opuFhat.cpp"
#include "opuUhat.cpp"
#include "opuStab.cpp"
#include "opuInitu.cpp"
#include "opuInitq.cpp"
#include "opuInitudg.cpp"
#include "opuInitwdg.cpp"
#include "opuInitodg.cpp"

#ifdef _ENZYME   
template <typename... Args>
void __enzyme_autodiff(void*, Args... args);
void __enzyme_fwddiff(void*, ...);
int enzyme_const, enzyme_dup;

template <typename T> void opuFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
        // printf("Calling opuFLUXEnzyme\n");
    __enzyme_fwddiff((void*)opuFlux<T>, 
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
template void opuFluxEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);

template <typename T> void opuSourceEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
    __enzyme_fwddiff((void*)opuSource<T>, 
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
template void opuSourceEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);


template <typename T> void opuUbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
    __enzyme_fwddiff((void*)opuUbou<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg, dudg,
                enzyme_const, odg,
                enzyme_dup, wdg, dwdg,
                enzyme_const, uhg,
                enzyme_const, nlg,
                enzyme_const, tau,
                enzyme_const, uinf,
                enzyme_const, param,
                enzyme_const, time,
                enzyme_const, modelnumber, 
                enzyme_const, ib,
                enzyme_const, ng,
                enzyme_const, nc,
                enzyme_const, ncu,
                enzyme_const, nd,
                enzyme_const, ncx,
                enzyme_const, nco,
                enzyme_const, ncw);            
}
template void opuUbouEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
        double, int, int, int, int, int, int, int, int, int);

template <typename T> void opuFbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *dodg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{    
    __enzyme_fwddiff((void*)opuFbou<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg, dudg,
                enzyme_dup, odg, dodg,
                enzyme_dup, wdg, dwdg,
                enzyme_dup, uhg, duhg,
                enzyme_const, nlg,
                enzyme_const, tau,
                enzyme_const, uinf,
                enzyme_const, param,
                enzyme_const, time,
                enzyme_const, modelnumber, 
                enzyme_const, ib,
                enzyme_const, ng,
                enzyme_const, nc,
                enzyme_const, ncu,
                enzyme_const, nd,
                enzyme_const, ncx,
                enzyme_const, nco,
                enzyme_const, ncw);            
}
template void opuFbouEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double  *, 
        double, int, int, int, int, int, int, int, int, int);

template <typename T> void opuFhatEnzyme(T *f, T *df, T *xg, T *udg1, T *dudg1, T *udg2, T *dudg2, 
        T *odg1, T *dodg1, T *odg2, T *dodg2, T *wdg1, T *dwdg1, T *wdg2, T *dwdg2, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
        __enzyme_fwddiff((void*)opuFhat<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg1, dudg1,
                enzyme_dup, udg2, dudg2,
                enzyme_dup, odg1, dodg1,
                enzyme_dup, odg2, dodg2,
                enzyme_dup, wdg1, dwdg1,
                enzyme_dup, wdg2, dwdg2,
                enzyme_dup, uhg, duhg, 
                enzyme_const, nlg,
                enzyme_const, tau,
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

template void opuFhatEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *,
        double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);


template <typename T> void opuStabEnzyme(T *f, T *df, T *xg, T *udg1, T *dudg1, T *udg2, T *dudg2,
        T *odg1, T *dodg1, T *odg2, T *dodg2, T *wdg1, T *dwdg1, T *wdg2, T *dwdg2, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param,
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
        __enzyme_fwddiff((void*)opuStab<T>, 
                enzyme_dup, f, df,
                enzyme_const, xg,
                enzyme_dup, udg1, dudg1,
                enzyme_dup, udg2, dudg2,
                enzyme_dup, odg1, dodg1,
                enzyme_dup, odg2, dodg2,
                enzyme_dup, wdg1, dwdg1,
                enzyme_dup, wdg2, dwdg2,
                enzyme_dup, uhg, duhg, 
                enzyme_const, nlg,
                enzyme_const, tau,
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
template void opuStabEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, double *,
        double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, 
        double, int, int, int, int, int, int, int, int);

template <typename T> void opuAvfieldEnzyme(T *f, T *df, T *xdg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, 
        T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne)
{
        __enzyme_fwddiff((void*)opuAvfield<T>, 
                enzyme_dup, f, df, 
                enzyme_const, xdg, 
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
template void opuAvfieldEnzyme(double *, double *, double *, double *, double *, double *, double *, double *, 
        double *, double *, double, int, int, int, int, int, int, int, int, int, int, int);
#endif

#endif

