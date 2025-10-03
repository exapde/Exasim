void KokkosUbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype odg1 = odg[0*ng+i];
		dstype odg2 = odg[1*ng+i];
		dstype odg3 = odg[2*ng+i];
		dstype odg4 = odg[3*ng+i];
		f[0*ng+i] = odg1;
		f[1*ng+i] = odg2;
		f[2*ng+i] = odg3;
		f[3*ng+i] = odg4;
	});
}

void KokkosUbou2(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou2", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
	});
}

void KokkosUbou3(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou3", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param9 = param[8];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype udg1 = udg[0*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = (param9*param11*udg1)/param10;
	});
}

void KokkosUbou4(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Ubou4", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype odg1 = odg[0*ng+i];
		dstype odg2 = odg[1*ng+i];
		dstype odg3 = odg[2*ng+i];
		dstype odg4 = odg[3*ng+i];
		f[0*ng+i] = odg1;
		f[1*ng+i] = odg2;
		f[2*ng+i] = odg3;
		f[3*ng+i] = odg4;
	});
}

void KokkosUbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		KokkosUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		KokkosUbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		KokkosUbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

