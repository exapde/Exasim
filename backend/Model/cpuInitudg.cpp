#include "cpuInitudg1.cpp"
#include "cpuInitudg2.cpp"

void cpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		cpuInitudg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitudg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

