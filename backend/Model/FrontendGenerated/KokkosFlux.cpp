void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype t2 = udg3+udg6;
		dstype t3 = udg4+udg5;
		dstype t4 = param1*t3;
		dstype t5 = param2*t2;
		f[0*ng+i] = t5+param1*udg3*2.0;
		f[1*ng+i] = t4;
		f[2*ng+i] = t4;
		f[3*ng+i] = t5+param1*udg6*2.0;
	});
}

