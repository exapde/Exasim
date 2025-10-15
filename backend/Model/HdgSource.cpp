void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		{
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = 0.0;
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
	});
}

