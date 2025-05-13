void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = 1.0/(udg1*udg1);
		dstype t4 = udg5*udg5;
		dstype t5 = t4*(1.0/2.0);
		dstype t6 = udg6*udg6;
		dstype t7 = t6*(1.0/2.0);
		dstype t8 = 1.0/udg1;
		dstype t9 = param1-1.0;
		dstype t10 = t2*t3*(1.0/2.0);
		dstype t11 = udg3*udg3;
		dstype t12 = t3*t11*(1.0/2.0);
		dstype t13 = t10+t12;
		dstype t14 = t13*udg1;
		dstype t15 = udg7*udg7;
		dstype t16 = t15*(1.0/2.0);
		dstype t17 = t5+t7+t14+t16-udg4;
		dstype t18 = t8*udg2*udg3;
		dstype t19 = t18-udg5*udg6;
		dstype t20 = t8*udg2*udg5;
		dstype t21 = t8*udg3*udg6;
		dstype t22 = t20+t21;
		dstype t23 = t8*udg4;
		dstype t24 = t5+t7;
		dstype t25 = t8*t24;
		dstype t26 = t23+t25-t8*t9*t17;
		dstype t27 = t8*udg2*udg6;
		f[0*ng+i] = udg2;
		f[1*ng+i] = -t5+t7+t2*t8-t9*t17;
		f[2*ng+i] = t19;
		f[3*ng+i] = -t22*udg5+t26*udg2+udg5*udg7;
		f[4*ng+i] = udg7;
		f[5*ng+i] = t27-t8*udg3*udg5;
		f[6*ng+i] = udg5;
		f[7*ng+i] = udg3;
		f[8*ng+i] = t19;
		f[9*ng+i] = t5-t7+t8*t11-t9*t17;
		f[10*ng+i] = -t22*udg6+t26*udg3+udg6*udg7;
		f[11*ng+i] = -t27+t8*udg3*udg5;
		f[12*ng+i] = udg7;
		f[13*ng+i] = udg6;
	});
}

