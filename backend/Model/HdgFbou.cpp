void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype uhg1 = uhg[0*ng+i];
		{
		f[0*ng+i] = -tau1*uhg1;
		}
		{
		f_udg[0*ng+i] = 0;
		f_udg[1*ng+i] = 0;
		f_udg[2*ng+i] = 0;
		}
		{
		f_uhg[0*ng+i] = -tau1;
		}
	});
}

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

