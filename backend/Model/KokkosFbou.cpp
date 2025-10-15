void KokkosFbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype tau1 = tau[0];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype t2 = udg2*udg2;
		dstype t3 = udg3*udg3;
		dstype t4 = param1-1.0;
		dstype t5 = 1.0/udg1;
		dstype t6 = t5*t5;
		dstype t7 = t5*udg4;
		dstype t8 = t2*t6;
		dstype t9 = t3*t6;
		dstype t10 = t8+t9;
		dstype t11 = (t10*udg1)/2.0;
		dstype t14 = -t4*(t11-udg4);
		dstype t15 = t5*t14;
		dstype t16 = t7+t15;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+tau1*(udg1-uhg1);
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*(t14+t2*t5)+nlg2*t5*udg2*udg3;
		f[2*ng+i] = tau1*(udg3-uhg3)+nlg2*(t14+t3*t5)+nlg1*t5*udg2*udg3;
		f[3*ng+i] = tau1*(udg4-uhg4)+nlg1*t16*udg2+nlg2*t16*udg3;
	});
}

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

