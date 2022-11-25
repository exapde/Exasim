#include "cpuInitudg1.cpp"
#include "cpuInitudg2.cpp"

template <typename T> void cpuInitudg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		cpuInitudg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitudg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void cpuInitudg(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitudg(float *, float *, float *, float *, int, int, int, int, int, int);
