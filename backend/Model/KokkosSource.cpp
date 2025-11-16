void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = udg3*udg3;
		dstype t3 = 1.0/3.141592653589793;
		dstype t4 = udg1*1.0E+2;
		dstype t5 = xdg2+1.0/1.0E+1;
		dstype t6 = udg1-1.0/2.0E+1;
		dstype t7 = t4-5.0;
		dstype t8 = 1.0/t5;
		dstype t9 = atan(t7);
		dstype t10 = t3*t9;
		dstype t11 = t10+1.0/2.0;
		dstype t12 = t6*t11;
		dstype t13 = t12+5.318299276490819E-2;
		dstype t14 = 1.0/t13;
		dstype t15 = t14*t14;
		f[0*ng+i] = -t8*udg3;
		f[1*ng+i] = -t8*t14*udg2*udg3;
		f[2*ng+i] = -t2*t8*t14;
		f[3*ng+i] = t8*udg3*((t2*t15)/5.0-t14*udg4*(7.0/5.0)+(t15*(udg2*udg2))/5.0);
	});
}

