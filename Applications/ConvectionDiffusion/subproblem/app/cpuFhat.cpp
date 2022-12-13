#include "cpuFhat1.cpp"
#include "cpuFhat2.cpp"

template <typename T> void cpuFhat(T *f, T *xdg, T *udg1, T *udg2,  T *odg1, T *odg2,  T *wdg1, T *wdg2,  T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (modelnumber == 1)
		cpuFhat1(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		cpuFhat2(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void cpuFhat(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double time, int, int, int, int, int, int, int, int);
template void cpuFhat(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float time, int, int, int, int, int, int, int, int);
