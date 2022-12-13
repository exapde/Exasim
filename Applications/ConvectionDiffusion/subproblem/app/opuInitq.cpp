#include "opuInitq1.cpp"
#include "opuInitq2.cpp"

template <typename T> void opuInitq(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		opuInitq1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		opuInitq2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void opuInitq(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitq(float *, float *, float *, float *, int, int, int, int, int, int);
