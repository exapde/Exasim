void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = param1-1.0;
		dstype t5 = 1.0/udg1;
		dstype t6 = t5*t5;
		dstype t7 = t5*udg4;
		dstype t8 = t5*udg2*udg3;
		dstype t9 = t2*t6;
		dstype t10 = t3*t6;
		dstype t11 = t9+t10;
		dstype t12 = (t11*udg1)/2.0;
		dstype t15 = -t4*(t12-udg4);
		dstype t16 = t5*t15;
		dstype t17 = t7+t16;
		f[0*ng+i] = udg2;
		f[1*ng+i] = t15+t2*t5;
		f[2*ng+i] = t8;
		f[3*ng+i] = t17*udg2;
		f[4*ng+i] = udg3;
		f[5*ng+i] = t8;
		f[6*ng+i] = t15+t3*t5;
		f[7*ng+i] = t17*udg3;
	});
}

