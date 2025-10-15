void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype uhg1 = uhg[0*ng+i];
		{
		dstype t2 = tau1*uhg1;
		dstype t3 = -t2;
		f[0*ng+i] = t3;
		f[1*ng+i] = t3;
		f[2*ng+i] = t3;
		f[3*ng+i] = t3;
		}
		{
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = t2;
		f_uhg[2*ng+i] = t2;
		f_uhg[3*ng+i] = t2;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = 0.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = 0.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = 0.0;
		}
	});
}

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

