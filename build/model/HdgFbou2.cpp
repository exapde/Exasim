void HdgFbou21(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou21", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		{
		f[0*ng+i] = -tau1*(uhg1-sin(xdg1*3.141592653589793)*sin(xdg2*3.141592653589793));
		}
		{
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -tau1;
		}
	});
}

void HdgFbou2(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou21(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

