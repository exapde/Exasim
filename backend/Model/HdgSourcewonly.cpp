#include "HdgSourcewonly1.cpp"
#include "HdgSourcewonly2.cpp"

void HdgSourcewonly(dstype* f, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (modelnumber == 1)
		HdgSourcewonly1(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		HdgSourcewonly2(f, f_wdg, xdg, udg, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

