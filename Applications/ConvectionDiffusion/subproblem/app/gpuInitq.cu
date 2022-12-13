#include "gpuInitq1.cu"
#include "gpuInitq2.cu"

template <typename T> void gpuInitq(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		gpuInitq1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		gpuInitq2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitq(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitq(float *, float *, float *, float *, int, int, int, int, int, int);
