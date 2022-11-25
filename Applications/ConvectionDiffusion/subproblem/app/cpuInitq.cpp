#include "cpuInitq1.cpp"
#include "cpuInitq2.cpp"

template <typename T> void cpuInitq(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		cpuInitq1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitq2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void cpuInitq(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitq(float *, float *, float *, float *, int, int, int, int, int, int);
