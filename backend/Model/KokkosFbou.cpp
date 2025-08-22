void KokkosFbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype param3 = param[2];
		dstype param4 = param[3];
		dstype param5 = param[4];
		dstype tau1 = tau[0];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype odg2 = odg[1*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = param5*param5;
		dstype t3 = xdg1*xdg1;
		dstype t4 = xdg2*xdg2;
		dstype t5 = 1.0/param2;
		dstype t6 = time-1.0;
		dstype t7 = -t2;
		dstype t8 = t3+t4+t7;
		dstype t9 = param4*t8;
		dstype t10 = cosh(t9);
		dstype t11 = 1.0/t10;
		dstype t12 = param3*t11;
		dstype t13 = t12+1.0;
		dstype t14 = t5*t6*t13;
		f[0*ng+i] = -nlg1*((odg1*udg1)/(t14-time)-param1*t12*udg2)-nlg2*((odg2*udg1)/(t14-time)-param1*t12*udg3)+tau1*(udg1-uhg1);
	});
}

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

