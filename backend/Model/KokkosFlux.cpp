void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype udg1 = udg[0*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype udg8 = udg[7*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype t2 = udg5+udg7;
		dstype t3 = udg6+udg10;
		dstype t4 = udg9+udg11;
		dstype t5 = udg4+udg8+udg12;
		dstype t6 = param1*t2;
		dstype t7 = param1*t3;
		dstype t8 = param1*t4;
		dstype t9 = param2*t5;
		f[0*ng+i] = t9+param1*udg4*2.0;
		f[1*ng+i] = t6;
		f[2*ng+i] = t7;
		f[3*ng+i] = t6;
		f[4*ng+i] = t9+param1*udg8*2.0;
		f[5*ng+i] = t8;
		f[6*ng+i] = t7;
		f[7*ng+i] = t8;
		f[8*ng+i] = t9+param1*udg12*2.0;
	});
}

