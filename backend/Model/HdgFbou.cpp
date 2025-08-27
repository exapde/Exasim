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
		dstype udg7 = udg[6*ng+i];
		dstype udg8 = udg[7*ng+i];
		dstype udg9 = udg[8*ng+i];
		dstype udg10 = udg[9*ng+i];
		dstype udg11 = udg[10*ng+i];
		dstype udg12 = udg[11*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype nlg3 = nlg[2*ng+i];
		{
		dstype t2 = udg5+udg7;
		dstype t3 = udg6+udg10;
		dstype t4 = udg9+udg11;
		dstype t5 = udg4+udg8+udg12;
		dstype t6 = param2*t5;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*(t6+param1*udg4*2.0)+nlg2*param1*t2+nlg3*param1*t3;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg2*(t6+param1*udg8*2.0)+nlg1*param1*t2+nlg3*param1*t4;
		f[2*ng+i] = tau1*(udg3-uhg3)+nlg3*(t6+param1*udg12*2.0)+nlg1*param1*t3+nlg2*param1*t4;
		}
		{
		dstype t2 = nlg1*param1;
		dstype t3 = nlg1*param2;
		dstype t4 = nlg2*param1;
		dstype t5 = nlg2*param2;
		dstype t6 = nlg3*param1;
		dstype t7 = nlg3*param2;
		dstype t8 = param1*2.0;
		dstype t9 = param2+t8;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = tau1;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = tau1;
		f_udg[9*ng+i] = nlg1*t9;
		f_udg[10*ng+i] = t5;
		f_udg[11*ng+i] = t7;
		f_udg[12*ng+i] = t4;
		f_udg[13*ng+i] = t2;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = t6;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = t2;
		f_udg[18*ng+i] = t4;
		f_udg[19*ng+i] = t2;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = t3;
		f_udg[22*ng+i] = nlg2*t9;
		f_udg[23*ng+i] = t7;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t6;
		f_udg[26*ng+i] = t4;
		f_udg[27*ng+i] = t6;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = t2;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = t6;
		f_udg[32*ng+i] = t4;
		f_udg[33*ng+i] = t3;
		f_udg[34*ng+i] = t5;
		f_udg[35*ng+i] = nlg3*t9;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = t2;
		f_uhg[5*ng+i] = 0.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = t2;
		}
	});
}

void HdgFbou2(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype tau1 = tau[0];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		{
		f[0*ng+i] = -tau1*uhg1;
		f[1*ng+i] = -tau1*uhg2;
		f[2*ng+i] = -tau1*uhg3;
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
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 0.0;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = 0.0;
		f_udg[18*ng+i] = 0.0;
		f_udg[19*ng+i] = 0.0;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = 0.0;
		f_udg[22*ng+i] = 0.0;
		f_udg[23*ng+i] = 0.0;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = 0.0;
		f_udg[26*ng+i] = 0.0;
		f_udg[27*ng+i] = 0.0;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = 0.0;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = 0.0;
		f_udg[32*ng+i] = 0.0;
		f_udg[33*ng+i] = 0.0;
		f_udg[34*ng+i] = 0.0;
		f_udg[35*ng+i] = 0.0;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = t2;
		f_uhg[5*ng+i] = 0.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = t2;
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
		dstype param5 = param[4];
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
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		dstype nlg3 = nlg[2*ng+i];
		{
		dstype t2 = udg5+udg7;
		dstype t3 = udg6+udg10;
		dstype t4 = udg9+udg11;
		dstype t5 = udg4+udg8+udg12;
		dstype t6 = param2*t5;
		f[0*ng+i] = -param3+tau1*(udg1-uhg1)+nlg1*(t6+param1*udg4*2.0)+nlg2*param1*t2+nlg3*param1*t3;
		f[1*ng+i] = -param4+tau1*(udg2-uhg2)+nlg2*(t6+param1*udg8*2.0)+nlg1*param1*t2+nlg3*param1*t4;
		f[2*ng+i] = -param5+tau1*(udg3-uhg3)+nlg3*(t6+param1*udg12*2.0)+nlg1*param1*t3+nlg2*param1*t4;
		}
		{
		dstype t2 = nlg1*param1;
		dstype t3 = nlg1*param2;
		dstype t4 = nlg2*param1;
		dstype t5 = nlg2*param2;
		dstype t6 = nlg3*param1;
		dstype t7 = nlg3*param2;
		dstype t8 = param1*2.0;
		dstype t9 = param2+t8;
		f_udg[0*ng+i] = tau1;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = tau1;
		f_udg[5*ng+i] = 0.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = tau1;
		f_udg[9*ng+i] = nlg1*t9;
		f_udg[10*ng+i] = t5;
		f_udg[11*ng+i] = t7;
		f_udg[12*ng+i] = t4;
		f_udg[13*ng+i] = t2;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = t6;
		f_udg[16*ng+i] = 0.0;
		f_udg[17*ng+i] = t2;
		f_udg[18*ng+i] = t4;
		f_udg[19*ng+i] = t2;
		f_udg[20*ng+i] = 0.0;
		f_udg[21*ng+i] = t3;
		f_udg[22*ng+i] = nlg2*t9;
		f_udg[23*ng+i] = t7;
		f_udg[24*ng+i] = 0.0;
		f_udg[25*ng+i] = t6;
		f_udg[26*ng+i] = t4;
		f_udg[27*ng+i] = t6;
		f_udg[28*ng+i] = 0.0;
		f_udg[29*ng+i] = t2;
		f_udg[30*ng+i] = 0.0;
		f_udg[31*ng+i] = t6;
		f_udg[32*ng+i] = t4;
		f_udg[33*ng+i] = t3;
		f_udg[34*ng+i] = t5;
		f_udg[35*ng+i] = nlg3*t9;
		}
		{
		dstype t2 = -tau1;
		f_uhg[0*ng+i] = t2;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = t2;
		f_uhg[5*ng+i] = 0.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = t2;
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

