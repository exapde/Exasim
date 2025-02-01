void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg1 = xdg[0*ng+i];
		dstype t2 = xdg1-1.0/2.0;
		f[0*ng+i] = (t2*t2)*-4.0+1.0;
	});
}

