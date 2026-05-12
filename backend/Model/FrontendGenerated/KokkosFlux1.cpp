void KokkosFlux1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		f[0*ng+i] = param1*udg2;
		f[1*ng+i] = param1*udg3;
	});
}

