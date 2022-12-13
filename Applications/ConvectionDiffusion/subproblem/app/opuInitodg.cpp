#include "opuInitodg1.cpp"
#include "opuInitodg2.cpp"

template <typename T> void opuInitodg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		opuInitodg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		opuInitodg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void opuInitodg(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitodg(float *, float *, float *, float *, int, int, int, int, int, int);
