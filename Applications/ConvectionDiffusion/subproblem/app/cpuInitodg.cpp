#include "cpuInitodg1.cpp"
#include "cpuInitodg2.cpp"

template <typename T> void cpuInitodg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		cpuInitodg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitodg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void cpuInitodg(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitodg(float *, float *, float *, float *, int, int, int, int, int, int);
