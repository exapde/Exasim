#include "HdgFint1.cpp"
#include "HdgFint2.cpp"

void HdgFint(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (modelnumber == 1)
		HdgFint1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		HdgFint2(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ib, ng, nc, ncu, nd, ncx, nco, ncw);
}

