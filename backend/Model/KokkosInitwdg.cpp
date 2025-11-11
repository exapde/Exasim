#include "KokkosInitwdg1.cpp"
#include "KokkosInitwdg2.cpp"

void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		KokkosInitwdg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		KokkosInitwdg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

