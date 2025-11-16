void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype udg8 = udg[7*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = 1.0/3.141592653589793;
		dstype t5 = udg1*1.0E+2;
		dstype t6 = udg1-1.0/2.0E+1;
		dstype t7 = t5-5.0;
		dstype t8 = atan(t7);
		dstype t9 = t4*t8;
		dstype t10 = t9+1.0/2.0;
		dstype t11 = t6*t10;
		dstype t12 = t11+5.318299276490819E-2;
		dstype t13 = 1.0/t12;
		dstype t14 = t13*t13;
		dstype t15 = t13*udg2*udg3;
		dstype t16 = t13*udg4*(2.0/5.0);
		dstype t17 = t13*udg4*(7.0/5.0);
		dstype t18 = -t16;
		dstype t19 = -t17;
		dstype t20 = (t2*t14)/5.0;
		dstype t21 = (t3*t14)/5.0;
		dstype t22 = t18+t20+t21;
		dstype t23 = t19+t20+t21;
		dstype t24 = t12*t22;
		dstype t25 = -t24;
		f[0*ng+i] = udg2+odg1*udg5;
		f[1*ng+i] = t25+odg1*udg6+t2*t13;
		f[2*ng+i] = t15+odg1*udg7;
		f[3*ng+i] = odg1*udg8-t23*udg2;
		f[4*ng+i] = udg3+odg1*udg9;
		f[5*ng+i] = t15+odg1*udg10;
		f[6*ng+i] = t25+odg1*udg11+t3*t13;
		f[7*ng+i] = odg1*udg12-t23*udg3;
	});
}

