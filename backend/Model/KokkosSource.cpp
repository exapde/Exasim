void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		f[0*ng+i] = sin(xdg1*3.141592653589793)*sin(xdg2*3.141592653589793)*1.973920880217872E+1;
	});
}

