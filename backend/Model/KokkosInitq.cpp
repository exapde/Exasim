#include "KokkosInitq1.cpp"
#include "KokkosInitq2.cpp"

void KokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		KokkosInitq1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		KokkosInitq2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

