void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param14 = param[13];
		dstype param16 = param[15];
		dstype param18 = param[17];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		{
		dstype t2 = 1.0/3.141592653589793;
		dstype t3 = udg2*1.0E+3;
		dstype t4 = atan(t3);
		dstype t5 = t2*t4;
		dstype t6 = t5+1.0/2.0;
		dstype t7 = t6*udg2;
		f[0*ng+i] = param18+tau1*(udg1-uhg1)+nlg1*udg3*(param16*(t7+3.183097800805168E-4)-param14*(t7-9.996816902199195E-1));
		f[1*ng+i] = -tau1*(uhg2-1.0);
		}
		{
		dstype t2 = udg2*udg2;
		dstype t3 = 1.0/3.141592653589793;
		dstype t4 = udg2*1.0E+3;
		dstype t5 = atan(t4);
		dstype t6 = t2*1.0E+6;
		dstype t7 = t3*t5;
		dstype t8 = t6+1.0;
		dstype t9 = t7+1.0/2.0;
		dstype t10 = 1.0/t8;
		dstype t11 = t9*udg2;
		dstype t12 = t3*t4*t10;
		dstype t13 = t9+t12;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = -nlg1*udg3*(param14*t13-param16*t13);
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = nlg1*(param16*(t11+3.183097800805168E-4)-param14*(t11-9.996816902199195E-1));
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
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
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param5 = param[4];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype param12 = param[11];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		{
		dstype t2 = 1.0/param3;
		f[0*ng+i] = -tau1*(uhg1-param2*t2);
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*1.0/(param1*param1)*param10*param12*udg4*exp(-(param11*t2)/(param5*udg1));
		}
		{
		dstype t2 = 1.0/(param1*param1);
		dstype t3 = 1.0/param3;
		dstype t4 = 1.0/param5;
		dstype t5 = 1.0/udg1;
		dstype t6 = param11*t3*t4*t5;
		dstype t7 = -t6;
		dstype t8 = exp(t7);
		f_udg[0*ng+i] = 0.0;
		f_udg[1*ng+i] = nlg1*param10*param12*t2*t5*t6*t8*udg4;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = tau1;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = nlg1*param10*param12*t2*t8;
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
}

