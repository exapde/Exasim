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
		dstype t2 = udg3+udg6;
		dstype t3 = udg4+udg5;
		dstype t4 = param2*t2;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*(t4+param1*udg3*2.0)+nlg2*param1*t3;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg2*(t4+param1*udg6*2.0)+nlg1*param1*t3;
		}
		{
		dstype t2 = nlg1*param1;
		dstype t3 = nlg2*param1;
		dstype t4 = param1*2.0;
		dstype t5 = param2+t4;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = tau1;
		f_udg[4*ng+i] = nlg1*t5;
		f_udg[5*ng+i] = nlg2*param2;
		f_udg[6*ng+i] = t3;
		f_udg[7*ng+i] = t2;
		f_udg[8*ng+i] = t3;
		f_udg[9*ng+i] = t2;
		f_udg[10*ng+i] = nlg1*param2;
		f_udg[11*ng+i] = nlg2*t5;
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
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		{
		f[0*ng+i] = -tau1*uhg1;
		f[1*ng+i] = -tau1*uhg2;
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
		dstype t2 = udg3+udg6;
		dstype t3 = udg4+udg5;
		dstype t4 = param2*t2;
		f[0*ng+i] = -param3+tau1*(udg1-uhg1)+nlg1*(t4+param1*udg3*2.0)+nlg2*param1*t3;
		f[1*ng+i] = -param4+tau1*(udg2-uhg2)+nlg2*(t4+param1*udg6*2.0)+nlg1*param1*t3;
		}
		{
		dstype t2 = nlg1*param1;
		dstype t3 = nlg2*param1;
		dstype t4 = param1*2.0;
		dstype t5 = param2+t4;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = tau1;
		f_udg[4*ng+i] = nlg1*t5;
		f_udg[5*ng+i] = nlg2*param2;
		f_udg[6*ng+i] = t3;
		f_udg[7*ng+i] = t2;
		f_udg[8*ng+i] = t3;
		f_udg[9*ng+i] = t2;
		f_udg[10*ng+i] = nlg1*param2;
		f_udg[11*ng+i] = nlg2*t5;
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

