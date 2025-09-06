void HdgFbou1(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		dstype odg2 = odg[1*ng+i];
		dstype odg3 = odg[2*ng+i];
		dstype odg4 = odg[3*ng+i];
		dstype odg5 = odg[4*ng+i];
		{
		f[0*ng+i] = odg2-uhg1;
		f[1*ng+i] = odg3-uhg2;
		f[2*ng+i] = odg4-uhg3;
		f[3*ng+i] = odg5-uhg4;
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
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou2(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
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
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou3(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou3", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = udg2-uhg2;
		f[2*ng+i] = udg3-uhg3;
		f[3*ng+i] = udg4-uhg4;
		}
		{
		f_udg[0*ng+i] = 1.0;
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = 1.0;
		f_udg[6*ng+i] = 0.0;
		f_udg[7*ng+i] = 0.0;
		f_udg[8*ng+i] = 0.0;
		f_udg[9*ng+i] = 0.0;
		f_udg[10*ng+i] = 1.0;
		f_udg[11*ng+i] = 0.0;
		f_udg[12*ng+i] = 0.0;
		f_udg[13*ng+i] = 0.0;
		f_udg[14*ng+i] = 0.0;
		f_udg[15*ng+i] = 1.0;
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
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = 0.0;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
		}
	});
}

void HdgFbou4(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou4", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		{
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = -uhg4+(param9*param11*uhg1)/param10;
		}
		{
		f_udg[0*ng+i] = 1.0;
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
		f_udg[36*ng+i] = 0.0;
		f_udg[37*ng+i] = 0.0;
		f_udg[38*ng+i] = 0.0;
		f_udg[39*ng+i] = 0.0;
		f_udg[40*ng+i] = 0.0;
		f_udg[41*ng+i] = 0.0;
		f_udg[42*ng+i] = 0.0;
		f_udg[43*ng+i] = 0.0;
		f_udg[44*ng+i] = 0.0;
		f_udg[45*ng+i] = 0.0;
		f_udg[46*ng+i] = 0.0;
		f_udg[47*ng+i] = 0.0;
		}
		{
		f_uhg[0*ng+i] = -1.0;
		f_uhg[1*ng+i] = 0.0;
		f_uhg[2*ng+i] = 0.0;
		f_uhg[3*ng+i] = (param9*param11)/param10;
		f_uhg[4*ng+i] = 0.0;
		f_uhg[5*ng+i] = -1.0;
		f_uhg[6*ng+i] = 0.0;
		f_uhg[7*ng+i] = 0.0;
		f_uhg[8*ng+i] = 0.0;
		f_uhg[9*ng+i] = 0.0;
		f_uhg[10*ng+i] = -1.0;
		f_uhg[11*ng+i] = 0.0;
		f_uhg[12*ng+i] = 0.0;
		f_uhg[13*ng+i] = 0.0;
		f_uhg[14*ng+i] = 0.0;
		f_uhg[15*ng+i] = -1.0;
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
	else if (ib == 4)
		HdgFbou4(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

