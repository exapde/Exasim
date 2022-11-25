#include "gpuInitudg1.cu"
#include "gpuInitudg2.cu"

template <typename T> void gpuInitudg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		gpuInitudg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		gpuInitudg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void gpuInitudg(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitudg(float *, float *, float *, float *, int, int, int, int, int, int);
