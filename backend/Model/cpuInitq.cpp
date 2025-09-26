#include "cpuInitq1.cpp"
#include "cpuInitq2.cpp"

void cpuInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		cpuInitq1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitq2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

