#include "opuInitu1.cpp"
#include "opuInitu2.cpp"

template <typename T> void opuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		opuInitu1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		opuInitu2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void opuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
