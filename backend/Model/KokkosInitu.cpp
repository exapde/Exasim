#include "KokkosInitu1.cpp"
#include "KokkosInitu2.cpp"

void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		KokkosInitu1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		KokkosInitu2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

