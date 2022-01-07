#ifndef __APPLICATION_H__
#define __APPLICATION_H__

#ifdef HAVE_ONETHREAD
template <typename T> void opuFlux(T *f, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuSource(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuSourcew(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void opuTdfunc(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuFbou(T *fh, T *pg, T *udg, T *odg, T *wdg,  T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void opuUbou(T *ub, T *pg, T *udg, T *odg, T *wdg, T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void opuAvfield(T *avfield, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void opuOutput(T *output, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void opuFhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuUhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuStab(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuInitu(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void opuInitq(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void opuInitudg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void opuInitwdg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int ncw, int npe, int ne);
template <typename T> void opuInitodg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nco, int npe, int ne);

#ifdef HAVE_ENZYME       
template <typename T> void opuFluxEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuSourceEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
// TODO: do the variable names need to be made consisten here...i.e. using nl instead of nlg
template <typename T> void opuUbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg,T *odg, T *wdg, T *dwdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuFbouEnzyme(T *f, T *df, T *xg, T *udg, T *dudg, T *odg, T *wdg, T *dwdg, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuFhatEnzyme(T *f, T *df, T *xg, T *udg1, T *dudg1, T *udg2, T *dudg2, 
        T *odg1, T *odg2,  T *wdg1, T *dwdg1, T *wdg2, T *dwdg2, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void opuStabEnzyme(T *f, T *df, T *xg, T *udg1, T *dudg1, T *udg2, T *dudg2, 
        T *odg1, T *odg2,  T *wdg1, T *dwdg1, T *wdg2, T *dwdg2, T *uhg, T *duhg, T *nlg, T *tau, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
#endif

#endif                

#ifdef HAVE_OPENMP        
template <typename T> void cpuFlux(T *f, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuSource(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuSourcew(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void cpuTdfunc(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuFbou(T *fh, T *pg, T *udg, T *odg, T *wdg,  T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void cpuUbou(T *ub, T *pg, T *udg, T *odg, T *wdg, T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void cpuAvfield(T *avfield, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void cpuOutput(T *output, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void cpuFhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuUhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuStab(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void cpuInitu(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void cpuInitq(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void cpuInitudg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void cpuInitwdg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int ncw, int npe, int ne);
template <typename T> void cpuInitodg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nco, int npe, int ne);
#endif                

#ifdef HAVE_CUDA      
template <typename T> void gpuFlux(T *f, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuSource(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuSourcew(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void gpuTdfunc(T *s, T *pg, T *udg, T *odg, T *wdg, T *uinf, T *param, 
        T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuFbou(T *fh, T *pg, T *udg, T *odg, T *wdg,  T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void gpuUbou(T *ub, T *pg, T *udg, T *odg, T *wdg, T *uhg, T *nl, 
        T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, 
        int ncx, int nco, int ncw);
template <typename T> void gpuAvfield(T *avfield, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void gpuOutput(T *output, T *xdg, T *udg, T *odg, T *wdg, T *uinf, 
        T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne);
template <typename T> void gpuFhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuUhat(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuStab(T *fh, T *xg, T *udg1, T *udg2, T *odg1, T *odg2,
        T *wdg1, T *wdg2, T *uhg, T *nl, T *tau, T *uinf, T *param, T time, int modelnumber, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw);
template <typename T> void gpuInitu(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void gpuInitq(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void gpuInitudg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nc, int npe, int ne);
template <typename T> void gpuInitwdg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int ncw, int npe, int ne);
template <typename T> void gpuInitodg(T *udg, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nco, int npe, int ne);
#endif                

#endif  

