#include "gpuFlux1.cu"
#include "gpuFlux2.cu"

template <typename T> void gpuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (modelnumber == 1)
		gpuFlux1(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		gpuFlux2(f, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void gpuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
