void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
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
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = param1-1.0;
		dstype t5 = 1.0/param2;
		dstype t6 = 1.0/param3;
		dstype t7 = 1.0/udg1;
		dstype t8 = t7*t7;
		dstype t9 = t7*udg4;
		dstype t10 = t7*udg2*udg3;
		dstype t11 = t7*udg2*udg5;
		dstype t12 = t7*udg3*udg5;
		dstype t13 = t7*udg2*udg9;
		dstype t14 = t7*udg3*udg9;
		dstype t15 = 1.0/t4;
		dstype t20 = (t2*t8)/2.0;
		dstype t21 = (t3*t8)/2.0;
		dstype t29 = -t7*(t14-udg11);
		dstype t31 = t7*(t14-udg11)*-2.0;
		dstype t43 = -t5*(t7*(t12-udg7)+t7*(t13-udg10));
		dstype t34 = t20+t21;
		dstype t45 = t10+t43;
		dstype t35 = t34*udg1;
		dstype t39 = -t4*(t35-udg4);
		dstype t42 = t7*t39;
		dstype t44 = t9+t42;
		f[0*ng+i] = udg2;
		f[1*ng+i] = t39-t5*(t29+t7*(t11-udg6)*2.0)*(2.0/3.0)+t2*t7;
		f[2*ng+i] = t45;
		f[3*ng+i] = t44*udg2+t7*t43*udg3-t5*t7*udg2*(t29+t7*(t11-udg6)*2.0)*(2.0/3.0)+param1*t5*t6*t8*t15*(t4*udg1*(udg8+udg1*(t8*udg2*(t11-udg6)+t8*udg3*(t12-udg7))-t34*udg5)+t4*udg5*(t35-udg4));
		f[4*ng+i] = udg3;
		f[5*ng+i] = t45;
		f[6*ng+i] = t39+t5*(t31+t7*(t11-udg6))*(2.0/3.0)+t3*t7;
		f[7*ng+i] = t44*udg3+t7*t43*udg2+t5*t7*udg3*(t31+t7*(t11-udg6))*(2.0/3.0)+param1*t5*t6*t8*t15*(t4*udg1*(udg12+udg1*(t8*udg2*(t13-udg10)+t8*udg3*(t14-udg11))-t34*udg9)+t4*udg9*(t35-udg4));
	});
}

