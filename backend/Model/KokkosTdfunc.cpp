void KokkosTdfunc(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Tdfunc", ng, KOKKOS_LAMBDA(const size_t i) {
		f[0*ng+i] = 1.0;
		f[1*ng+i] = 1.0;
		f[2*ng+i] = 1.0;
		f[3*ng+i] = 1.0;
		f[4*ng+i] = 1.0;
	});
}

