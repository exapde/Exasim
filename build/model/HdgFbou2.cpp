void HdgFbou21(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou21", ng, KOKKOS_LAMBDA(const size_t i) {
<<<<<<< Updated upstream
		dstype tau1 = tau[0];
		dstype uhg1 = uhg[0*ng+i];
		f[0*ng+i] = -tau1*uhg1;
=======
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
	});
}

void HdgFbou22(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou22", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
	});
}

void HdgFbou23(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou23", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = -uhg4+(param9*param11*uhg1)/param10;
	});
}

void HdgFbou24(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou24", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = -uhg2;
		f[2*ng+i] = -uhg3;
		f[3*ng+i] = -uhg4+(param9*param11*uhg1)/param10;
	});
}

void HdgFbou25(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou25", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param5 = param[4];
		dstype param6 = param[5];
		dstype param7 = param[6];
		dstype param8 = param[7];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = param5-uhg1;
		f[1*ng+i] = param6-uhg2;
		f[2*ng+i] = param7-uhg3;
		f[3*ng+i] = param8-uhg4;
	});
}

void HdgFbou26(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou26", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype uhg1 = uhg[0*ng+i];
		dstype uhg2 = uhg[1*ng+i];
		dstype uhg3 = uhg[2*ng+i];
		dstype uhg4 = uhg[3*ng+i];
		f[0*ng+i] = udg1-uhg1;
		f[1*ng+i] = udg2-uhg2;
		f[2*ng+i] = udg3-uhg3;
		f[3*ng+i] = udg4-uhg4;
>>>>>>> Stashed changes
	});
}

void HdgFbou2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		HdgFbou21(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
<<<<<<< Updated upstream
=======
	else if (ib == 2)
		HdgFbou22(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		HdgFbou23(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		HdgFbou24(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		HdgFbou25(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 6)
		HdgFbou26(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
>>>>>>> Stashed changes
}

