#include "opuInitudg1.cpp"
#include "opuInitudg2.cpp"

template <typename T> void opuInitudg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		opuInitudg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		opuInitudg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void opuInitudg(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitudg(float *, float *, float *, float *, int, int, int, int, int, int);
