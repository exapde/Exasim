void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		{
		f[0*ng+i] = 2*pow(M_PI, 2)*sin(M_PI*xdg1)*sin(M_PI*xdg2);
		}
		{
		f_udg[0*ng+i] = 0;
		f_udg[1*ng+i] = 0;
		f_udg[2*ng+i] = 0;
		}
	});
}

