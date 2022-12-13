#include "cpuInitwdg1.cpp"
#include "cpuInitwdg2.cpp"

template <typename T> void cpuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		cpuInitwdg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitwdg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void cpuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int);
