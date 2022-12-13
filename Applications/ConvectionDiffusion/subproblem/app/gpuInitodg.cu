#include "gpuInitodg1.cu"
#include "gpuInitodg2.cu"

template <typename T> void gpuInitodg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		gpuInitodg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		gpuInitodg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitodg(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitodg(float *, float *, float *, float *, int, int, int, int, int, int);
