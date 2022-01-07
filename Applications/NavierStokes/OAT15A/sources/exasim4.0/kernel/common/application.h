#ifndef __APPLICATION_H__
#define __APPLICATION_H__

#ifdef HAVE_ONETHREAD
template <typename T> void opuFlux(T *f, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuSource(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuTdfunc(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuFbou(T *fh, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuUbou(T *ub, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuAVfield(T *avField, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne);
template <typename T> void opuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nl, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void opuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt, int ib, int jb, 
        int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);

#endif                

#ifdef HAVE_OPENMP        
template <typename T> void cpuFlux(T *f, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuSource(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuTdfunc(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuFbou(T *fh, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuUbou(T *ub, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuAVfield(T *avField, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne);
template <typename T> void cpuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nl, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void cpuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt, int ib, int jb, 
        int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
#endif                

#ifdef HAVE_CUDA      
template <typename T> void gpuFlux(T *f, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuSource(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuTdfunc(T *s, T *pg, T *udg, T *odg, T *param, T time, 
            int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuFbou(T *fh, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuUbou(T *ub, T *pg, T *udg, T *uhg, T *odg, T *nl, T *tau, T *uinf, T *param, 
        T time, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuAVfield(T *avField, T *xdg, T *udg, T *odg, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int npe, int ne);
template <typename T> void gpuStab(T *fh, T *xg, T *uhg, T *udg1, T *udg2, T *odg1, T *odg2, T *nl, T *param, T time, 
        int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
template <typename T> void gpuUboutdep(T *ub, T *xg, T *udg, T *uhg, T *odg, T *udg0, T *udg1, T *udg2, 
        T *uhg0, T *uhg1, T *uhg2, T *nl, T *tau, T *uinf, T *param, T time, T dt, int ib, int jb, 
        int tstage, int ng, int nc, int ncu, int nd, int ncx, int nco, int appname);
#endif                

#endif  

