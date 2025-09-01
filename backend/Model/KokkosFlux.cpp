void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		f[0*ng+i] = param1*udg1*xdg2;
		f[1*ng+i] = -param2*udg1*xdg1;
	});
}

