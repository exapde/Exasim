void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param5 = param[4];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype odg2 = odg[1*ng+i];
		{
		dstype t2 = param5*param5;
		dstype t3 = xdg1*xdg1;
		dstype t4 = xdg2*xdg2;
		dstype t5 = 1.0/param2;
		dstype t6 = time-1.0;
		dstype t7 = -t2;
		dstype t8 = t3+t4+t7;
		dstype t9 = param4*t8;
		dstype t10 = cosh(t9);
		dstype t11 = sinh(t9);
		dstype t12 = 1.0/t10;
		dstype t13 = t12*t12;
		dstype t14 = param3*t12;
		dstype t15 = t14+1.0;
		dstype t16 = t5*t6*t15;
		dstype t19 = -1.0/(t16-time);
		f[0*ng+i] = t19*udg1*(t5*t15+t5*t6*((odg1*param3*param4*t11*t13*xdg1*2.0)/(t16-time)+(odg2*param3*param4*t11*t13*xdg2*2.0)/(t16-time))-1.0);
		}
		{
		dstype t2 = param5*param5;
		dstype t3 = xdg1*xdg1;
		dstype t4 = xdg2*xdg2;
		dstype t5 = 1.0/param2;
		dstype t6 = time-1.0;
		dstype t7 = -t2;
		dstype t8 = t3+t4+t7;
		dstype t9 = param4*t8;
		dstype t10 = cosh(t9);
		dstype t11 = sinh(t9);
		dstype t12 = 1.0/t10;
		dstype t13 = t12*t12;
		dstype t14 = param3*t12;
		dstype t15 = t14+1.0;
		dstype t16 = t5*t6*t15;
		dstype t19 = -1.0/(t16-time);
		f_udg[0*ng+i] = t19*(t5*t15+t5*t6*((odg1*param3*param4*t11*t13*xdg1*2.0)/(t16-time)+(odg2*param3*param4*t11*t13*xdg2*2.0)/(t16-time))-1.0);
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		}
	});
}

