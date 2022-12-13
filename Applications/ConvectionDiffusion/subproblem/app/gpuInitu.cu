#include "gpuInitu1.cu"
#include "gpuInitu2.cu"

template <typename T> void gpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		gpuInitu1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		gpuInitu2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
