#include "KokkosInitudg1.cpp"
#include "KokkosInitudg2.cpp"

void KokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		KokkosInitudg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		KokkosInitudg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

