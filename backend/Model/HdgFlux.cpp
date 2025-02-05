void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		{
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
		}
		{
		dstype t2 = 1.0/(udg1*udg1);
		dstype t3 = udg2*udg2;
		dstype t4 = 1.0/(udg1*udg1*udg1);
		dstype t5 = udg3*udg3;
		dstype t6 = param1-1.0;
		dstype t7 = t2*t3*(1.0/2.0);
		dstype t8 = t2*t5*(1.0/2.0);
		dstype t9 = t3*t4*2.0;
		dstype t10 = t4*t5*2.0;
		dstype t11 = t9+t10;
		dstype t14 = t11*udg1*(1.0/2.0);
		dstype t12 = t7+t8-t14;
		dstype t13 = t2*t5;
		dstype t15 = t2*udg4;
		dstype t16 = t2*t3;
		dstype t17 = t13+t16;
		dstype t23 = t17*udg1*(1.0/2.0);
		dstype t18 = -t23+udg4;
		dstype t19 = t2*t6*t18;
		dstype t20 = 1.0/udg1;
		dstype t21 = t6*t12*t20;
		dstype t22 = t15+t19+t21;
		dstype t24 = t20*udg3;
		dstype t25 = t20*udg2;
		dstype t26 = t20*udg4;
		dstype t27 = t6*t18*t20;
		dstype t28 = t6*t20;
		dstype t29 = t20+t28;
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = -t2*t3-t6*t12;
		f_udg[2*ng+i] = -t2*udg2*udg3;
		f_udg[3*ng+i] = -t22*udg2;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = -t2*udg2*udg3;
		f_udg[6*ng+i] = -t13-t6*t12;
		f_udg[7*ng+i] = -t22*udg3;
		f_udg[8*ng+i] = 1.0;
		f_udg[9*ng+i] = t20*udg2*2.0-t6*t20*udg2;
		f_udg[10*ng+i] = t24;
		f_udg[11*ng+i] = t26+t27-t2*t3*t6;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = t24;
		f_udg[14*ng+i] = -t6*t20*udg2;
		f_udg[15*ng+i] = -t2*t6*udg2*udg3;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = -t6*t20*udg3;
		f_udg[18*ng+i] = t25;
		f_udg[19*ng+i] = -t2*t6*udg2*udg3;
		f_udg[20*ng+i] = 1.0;
		f_udg[21*ng+i] = t25;
		f_udg[22*ng+i] = t20*udg3*2.0-t6*t20*udg3;
		f_udg[23*ng+i] = t26+t27-t2*t5*t6;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t6;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = t29*udg2;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = t6;
		f_udg[31*ng+i] = t29*udg3;
		}
		{
		f_wdg[0*ng+i] = 0.0;
		f_wdg[1*ng+i] = 0.0;
		f_wdg[2*ng+i] = 0.0;
		f_wdg[3*ng+i] = 0.0;
		f_wdg[4*ng+i] = 0.0;
		f_wdg[5*ng+i] = 0.0;
		f_wdg[6*ng+i] = 0.0;
		f_wdg[7*ng+i] = 0.0;
		}
	});
}

