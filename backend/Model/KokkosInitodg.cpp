#include "KokkosInitodg1.cpp"
#include "KokkosInitodg2.cpp"

void KokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		KokkosInitodg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		KokkosInitodg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

