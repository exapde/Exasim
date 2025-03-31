#include "KokkosFhat1.cpp"
#include "KokkosFhat2.cpp"

void KokkosFhat(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (modelnumber == 1)
		KokkosFhat1(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (modelnumber == 2)
		KokkosFhat2(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

