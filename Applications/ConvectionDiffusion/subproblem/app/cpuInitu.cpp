#include "cpuInitu1.cpp"
#include "cpuInitu2.cpp"

template <typename T> void cpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		cpuInitu1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitu2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void cpuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
