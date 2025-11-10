void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = 1.0/udg1;
		dstype t5 = t4*t4;
		dstype t6 = t4*udg2*udg3;
		dstype t7 = t4*udg4*(2.0/5.0);
		dstype t8 = t4*udg4*(7.0/5.0);
		dstype t9 = -t7;
		dstype t10 = -t8;
		dstype t11 = (t2*t5)/5.0;
		dstype t12 = (t3*t5)/5.0;
		dstype t13 = t9+t11+t12;
		dstype t14 = t10+t11+t12;
		dstype t15 = t13*udg1;
		dstype t16 = -t15;
		f[0*ng+i] = udg2;
		f[1*ng+i] = t16+t2*t4;
		f[2*ng+i] = t6;
		f[3*ng+i] = -t14*udg2;
		f[4*ng+i] = udg3;
		f[5*ng+i] = t6;
		f[6*ng+i] = t16+t3*t4;
		f[7*ng+i] = -t14*udg3;
	});
}

