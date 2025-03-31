#include "cpuInitwdg1.cpp"
#include "cpuInitwdg2.cpp"

void cpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		cpuInitwdg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitwdg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

