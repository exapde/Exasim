#include "cpuInitu1.cpp"
#include "cpuInitu2.cpp"

void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		cpuInitu1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitu2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

