void KokkosFbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/udg1;
		dstype t3 = udg2*udg2;
		dstype t4 = 1.0/(udg1*udg1);
		dstype t5 = 1.0/param2;
		dstype t6 = udg7-t2*udg3*udg5;
		dstype t7 = t2*t6;
		dstype t8 = udg10-t2*udg2*udg9;
		dstype t9 = t2*t8;
		dstype t10 = t7+t9;
		dstype t11 = t5*t10;
		dstype t12 = t2*udg2*udg3;
		dstype t13 = t11+t12;
		dstype t14 = udg3*udg3;
		dstype t15 = t3*t4*(1.0/2.0);
		dstype t16 = t4*t14*(1.0/2.0);
		dstype t17 = t15+t16;
		dstype t18 = udg4-t17*udg1;
		dstype t19 = param1-1.0;
		dstype t20 = t18*t19;
		dstype t21 = udg6-t2*udg2*udg5;
		dstype t22 = udg11-t2*udg3*udg9;
		f[0*ng+i] = 0.0;
		f[1*ng+i] = nlg2*t13+tau1*(udg2-uhg2)+nlg1*(t20+t2*t3+t5*(t2*t21*2.0-t2*t22)*(2.0/3.0));
		f[2*ng+i] = nlg1*t13+tau1*(udg3-uhg3)+nlg2*(t20+t2*t14-t5*(t2*t21-t2*t22*2.0)*(2.0/3.0));
		f[3*ng+i] = 0.0;
	});
}

void KokkosFbou2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype tau1 = tau[0];
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
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/udg1;
		dstype t3 = udg2*udg2;
		dstype t4 = 1.0/(udg1*udg1);
		dstype t5 = 1.0/param2;
		dstype t24 = t2*udg3*udg5;
		dstype t6 = -t24+udg7;
		dstype t7 = t2*t6;
		dstype t25 = t2*udg2*udg9;
		dstype t8 = -t25+udg10;
		dstype t9 = t2*t8;
		dstype t10 = t7+t9;
		dstype t11 = t5*t10;
		dstype t12 = t2*udg2*udg3;
		dstype t13 = t11+t12;
		dstype t14 = udg3*udg3;
		dstype t15 = t3*t4*(1.0/2.0);
		dstype t16 = t4*t14*(1.0/2.0);
		dstype t17 = t15+t16;
		dstype t23 = t17*udg1;
		dstype t18 = -t23+udg4;
		dstype t19 = param1-1.0;
		dstype t20 = t18*t19;
		dstype t26 = t2*udg2*udg5;
		dstype t21 = -t26+udg6;
		dstype t28 = t2*udg3*udg9;
		dstype t22 = -t28+udg11;
		dstype t27 = t2*t21*2.0;
		dstype t29 = t27-t2*t22;
		dstype t30 = t2*udg4;
		dstype t31 = t2*t18*t19;
		dstype t32 = t30+t31;
		dstype t33 = t2*t21;
		dstype t34 = t33-t2*t22*2.0;
		dstype t35 = 1.0/param3;
		dstype t36 = 1.0/t19;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+tau1*(udg1-uhg1);
		f[1*ng+i] = nlg1*(t20+t2*t3+t5*t29*(2.0/3.0))+nlg2*t13+tau1*(udg2-uhg2);
		f[2*ng+i] = nlg2*(t20+t2*t14-t5*t34*(2.0/3.0))+nlg1*t13+tau1*(udg3-uhg3);
		f[3*ng+i] = nlg1*(t32*udg2+t2*t5*t10*udg3+t2*t5*t29*udg2*(2.0/3.0)-param1*t4*t5*t35*t36*(t19*udg1*(-udg8+t17*udg5+udg1*(t4*t6*udg3+t4*t21*udg2))+t18*t19*udg5))+nlg2*(t32*udg3+t2*t5*t10*udg2-t2*t5*t34*udg3*(2.0/3.0)-param1*t4*t5*t35*t36*(t19*udg1*(-udg12+t17*udg9+udg1*(t4*t8*udg2+t4*t22*udg3))+t18*t19*udg9))+tau1*(udg4-uhg4);
	});
}

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		KokkosFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

