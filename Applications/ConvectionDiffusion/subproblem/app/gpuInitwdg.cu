#include "gpuInitwdg1.cu"
#include "gpuInitwdg2.cu"

template <typename T> void gpuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		gpuInitwdg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		gpuInitwdg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int);
