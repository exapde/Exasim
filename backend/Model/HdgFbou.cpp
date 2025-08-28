void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
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
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = -t2;
		dstype t5 = -t3;
		dstype t6 = t2+t5;
		dstype t7 = t3+t4+1.0;
		dstype t8 = 1.0/t6;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg2*(param1*udg5+param1*t8*udg4+param2*t7*udg4)-nlg1*(-param1*udg3+param1*t8*udg6+param2*t7*udg6);
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*(param1*udg4+param1*t8*udg5+param2*t7*udg5)-nlg2*(-param1*udg6+param1*t8*udg3+param2*t7*udg3);
		}
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = udg3*udg3;
		dstype t5 = udg4*udg4;
		dstype t6 = udg5*udg5;
		dstype t7 = udg6*udg6;
		dstype t8 = param2*udg3*udg4;
		dstype t9 = param2*udg3*udg5;
		dstype t12 = param2*udg4*udg6;
		dstype t13 = param2*udg5*udg6;
		dstype t10 = param2*t2;
		dstype t11 = param2*t3;
		dstype t14 = -t2;
		dstype t15 = -t3;
		dstype t17 = t2+t15;
		dstype t18 = t3+t14+1.0;
		dstype t19 = 1.0/t17;
		dstype t21 = param2*t18;
		dstype t20 = t19*t19;
		dstype t22 = param1*t19;
		dstype t23 = param1*t20*udg3*udg4;
		dstype t24 = param1*t20*udg3*udg5;
		dstype t25 = param1*t2*t20;
		dstype t26 = param1*t3*t20;
		dstype t27 = param1*t20*udg4*udg6;
		dstype t28 = param1*t20*udg5*udg6;
		dstype t30 = t8+t23;
		dstype t31 = t9+t24;
		dstype t32 = t12+t27;
		dstype t33 = t13+t28;
		dstype t38 = t11+t21+t22+t26;
		dstype t34 = nlg2*t30;
		dstype t35 = nlg1*t33;
		dstype t36 = -t34;
		dstype t37 = -t35;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = tau1;
		f_udg[4*ng+i] = -nlg2*t32+nlg1*(param1+param2*t7+param1*t7*t20);
		f_udg[5*ng+i] = t37+nlg2*(t10-t21-t22+t25);
		f_udg[6*ng+i] = t37+nlg2*t38;
		f_udg[7*ng+i] = -nlg2*t31+nlg1*(param1+param2*t6+param1*t6*t20);
		f_udg[8*ng+i] = -nlg1*t32+nlg2*(param1+param2*t5+param1*t5*t20);
		f_udg[9*ng+i] = t36+nlg1*t38;
		f_udg[10*ng+i] = t36+nlg1*(t10-t21-t22+t25);
		f_udg[11*ng+i] = -nlg1*t31+nlg2*(param1+param2*t4+param1*t4*t20);
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = t2;
		}
	});
}

void HdgFbou2(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		{
		f[0*ng+i] = -tau1*(uhg1-xdg1);
		f[1*ng+i] = -tau1*(uhg2-xdg2);
		}
		{
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 0.0;
		f_udg[11*ng+i] = 0.0;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = t2;
		}
	});
}

void HdgFbou3(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou3", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype udg6 = udg[5*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = -t2;
		dstype t5 = -t3;
		dstype t6 = t2+t5;
		dstype t7 = t3+t4+1.0;
		dstype t8 = 1.0/t6;
		f[0*ng+i] = param3+tau1*(udg1-uhg1)+nlg2*(param1*udg5+param1*t8*udg4+param2*t7*udg4)-nlg1*(-param1*udg3+param1*t8*udg6+param2*t7*udg6);
		f[1*ng+i] = param4+tau1*(udg2-uhg2)+nlg1*(param1*udg4+param1*t8*udg5+param2*t7*udg5)-nlg2*(-param1*udg6+param1*t8*udg3+param2*t7*udg3);
		}
		{
		dstype t2 = udg3*udg6;
		dstype t3 = udg4*udg5;
		dstype t4 = udg3*udg3;
		dstype t5 = udg4*udg4;
		dstype t6 = udg5*udg5;
		dstype t7 = udg6*udg6;
		dstype t8 = param2*udg3*udg4;
		dstype t9 = param2*udg3*udg5;
		dstype t12 = param2*udg4*udg6;
		dstype t13 = param2*udg5*udg6;
		dstype t10 = param2*t2;
		dstype t11 = param2*t3;
		dstype t14 = -t2;
		dstype t15 = -t3;
		dstype t17 = t2+t15;
		dstype t18 = t3+t14+1.0;
		dstype t19 = 1.0/t17;
		dstype t21 = param2*t18;
		dstype t20 = t19*t19;
		dstype t22 = param1*t19;
		dstype t23 = param1*t20*udg3*udg4;
		dstype t24 = param1*t20*udg3*udg5;
		dstype t25 = param1*t2*t20;
		dstype t26 = param1*t3*t20;
		dstype t27 = param1*t20*udg4*udg6;
		dstype t28 = param1*t20*udg5*udg6;
		dstype t30 = t8+t23;
		dstype t31 = t9+t24;
		dstype t32 = t12+t27;
		dstype t33 = t13+t28;
		dstype t38 = t11+t21+t22+t26;
		dstype t34 = nlg2*t30;
		dstype t35 = nlg1*t33;
		dstype t36 = -t34;
		dstype t37 = -t35;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = tau1;
		f_udg[4*ng+i] = -nlg2*t32+nlg1*(param1+param2*t7+param1*t7*t20);
		f_udg[5*ng+i] = t37+nlg2*(t10-t21-t22+t25);
		f_udg[6*ng+i] = t37+nlg2*t38;
		f_udg[7*ng+i] = -nlg2*t31+nlg1*(param1+param2*t6+param1*t6*t20);
		f_udg[8*ng+i] = -nlg1*t32+nlg2*(param1+param2*t5+param1*t5*t20);
		f_udg[9*ng+i] = t36+nlg1*t38;
		f_udg[10*ng+i] = t36+nlg1*(t10-t21-t22+t25);
		f_udg[11*ng+i] = -nlg1*t31+nlg2*(param1+param2*t4+param1*t4*t20);
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = t2;
		}
	});
}

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		HdgFbou2(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		HdgFbou3(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

