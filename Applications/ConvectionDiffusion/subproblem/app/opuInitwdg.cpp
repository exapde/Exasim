#include "opuInitwdg1.cpp"
#include "opuInitwdg2.cpp"

template <typename T> void opuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	if (modelnumber == 1)
		opuInitwdg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		opuInitwdg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

template void opuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int);
template void opuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int);
