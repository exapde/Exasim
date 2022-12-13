#include "opuFlux1.cpp"
#include "opuFlux2.cpp"

template <typename T> void opuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (modelnumber == 1)
		opuFlux1(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		opuFlux2(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
