void KokkosFbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype udg7 = udg[6*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype uhg5 = uhg[4*ng+i];
		dstype uhg6 = uhg[5*ng+i];
		dstype uhg7 = uhg[6*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = 1.0/udg1;
		dstype t3 = udg2*udg2;
		dstype t4 = 1.0/(udg1*udg1);
		dstype t5 = udg5*udg5;
		dstype t6 = t5*(1.0/2.0);
		dstype t7 = udg6*udg6;
		dstype t8 = t7*(1.0/2.0);
		dstype t9 = udg5*udg6;
		dstype t10 = t9-t2*udg2*udg3;
		dstype t11 = udg3*udg3;
		dstype t12 = param1-1.0;
		dstype t13 = t3*t4*(1.0/2.0);
		dstype t14 = t4*t11*(1.0/2.0);
		dstype t15 = t13+t14;
		dstype t16 = t15*udg1;
		dstype t17 = udg7*udg7;
		dstype t18 = t17*(1.0/2.0);
		dstype t19 = t6+t8+t16+t18-udg4;
		dstype t20 = t12*t19;
		dstype t21 = t2*udg2*udg5;
		dstype t22 = t2*udg3*udg6;
		dstype t23 = t21+t22;
		dstype t24 = t2*udg4;
		dstype t25 = t6+t8;
		dstype t26 = t2*t25;
		dstype t27 = t24+t26-t2*t12*t19;
		dstype t28 = t2*udg2*udg6;
		dstype t29 = t28-t2*udg3*udg5;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+tau1*(udg1-uhg1);
		f[1*ng+i] = -nlg1*(t6-t8+t20-t2*t3)-nlg2*t10+tau1*(udg2-uhg2);
		f[2*ng+i] = -nlg1*t10+tau1*(udg3-uhg3)+nlg2*(t6-t8-t20+t2*t11);
		f[3*ng+i] = nlg1*(-t23*udg5+t27*udg2+udg5*udg7)+nlg2*(-t23*udg6+t27*udg3+udg6*udg7)+tau1*(udg4-uhg4);
		f[4*ng+i] = -nlg2*t29+nlg1*udg7+tau1*(udg5-uhg5);
		f[5*ng+i] = nlg1*t29+nlg2*udg7+tau1*(udg6-uhg6);
		f[6*ng+i] = nlg1*udg5+nlg2*udg6+tau1*(udg7-uhg7);
	});
}

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

