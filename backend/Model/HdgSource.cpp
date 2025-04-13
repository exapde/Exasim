void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype param9 = param[8];
		dstype param12 = param[11];
		dstype param15 = param[14];
		dstype param17 = param[16];
		dstype param19 = param[18];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		{
		dstype t2 = 1.0/3.141592653589793;
		dstype t3 = udg2*1.0E+3;
		dstype t4 = atan(t3);
		dstype t5 = t2*t4;
		dstype t6 = t5+1.0/2.0;
		dstype t7 = t6*udg2;
		dstype t8 = t7+3.183097800805168E-4;
		dstype t10 = t7-9.996816902199195E-1;
		dstype t9 = param17*t8;
		dstype t11 = param15*t10;
		dstype t12 = -t11;
		dstype t13 = t9+t12;
		dstype t14 = 1.0/t13;
		f[0*ng+i] = (param12*param19*t8*(param6*param7*t14+param6*param9*t14*(param8+param3*param4*udg1*log(t8))))/param3;
		f[1*ng+i] = -param6*param12*t8;
		}
		{
		dstype t2 = udg2*udg2;
		dstype t3 = 1.0/3.141592653589793;
		dstype t4 = 1.0/param3;
		dstype t5 = udg2*1.0E+3;
		dstype t6 = atan(t5);
		dstype t7 = t2*1.0E+6;
		dstype t8 = t3*t6;
		dstype t9 = t7+1.0;
		dstype t10 = t8+1.0/2.0;
		dstype t11 = 1.0/t9;
		dstype t12 = t10*udg2;
		dstype t13 = t3*t5*t11;
		dstype t14 = t10+t13;
		dstype t18 = t12+3.183097800805168E-4;
		dstype t21 = t12-9.996816902199195E-1;
		dstype t15 = param15*t14;
		dstype t16 = param17*t14;
		dstype t19 = log(t18);
		dstype t20 = param17*t18;
		dstype t23 = param15*t21;
		dstype t17 = -t16;
		dstype t22 = param3*param4*t19*udg1;
		dstype t25 = -t23;
		dstype t24 = param8+t22;
		dstype t26 = t15+t17;
		dstype t27 = t20+t25;
		dstype t28 = 1.0/t27;
		dstype t29 = t28*t28;
		f_udg[0*ng+i] = param4*param6*param9*param12*param19*t18*t19*t28;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = param12*param19*t4*t14*(param6*param7*t28+param6*param9*t24*t28)+param12*param19*t4*t18*(param6*param7*t26*t29+param6*param9*t24*t26*t29+(param3*param4*param6*param9*t14*t28*udg1)/t18);
		f_udg[3*ng+i] = -param6*param12*t14;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		}
	});
}

