void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = 1.0/(udg1*udg1);
		dstype t4 = 1.0/udg1;
		dstype t5 = t2*t3;
		dstype t6 = udg3*udg3;
		dstype t7 = t3*t6;
		dstype t8 = t5+t7;
		dstype t12 = t8*udg1*(1.0/2.0);
		dstype t9 = -t12+udg4;
		dstype t10 = param1-1.0;
		dstype t11 = t4*udg2*udg3;
		dstype t13 = t9*t10;
		dstype t14 = t4*udg4;
		dstype t15 = t4*t9*t10;
		dstype t16 = t14+t15;
		f[0*ng+i] = udg2;
		f[1*ng+i] = t13+t2*t4;
		f[2*ng+i] = t11;
		f[3*ng+i] = t16*udg2;
		f[4*ng+i] = udg3;
		f[5*ng+i] = t11;
		f[6*ng+i] = t13+t4*t6;
		f[7*ng+i] = t16*udg3;
	});
}

