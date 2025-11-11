#include "cpuInitodg1.cpp"
#include "cpuInitodg2.cpp"

void cpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	if (modelnumber == 1)
		cpuInitodg1(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
	else if (modelnumber == 2)
		cpuInitodg2(f, xdg, uinf, param, modelnumber, ng, ncx, nce, npe, ne);
}

