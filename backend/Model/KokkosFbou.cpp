void KokkosFbou1(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Fbou1", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype tau1 = tau[0];
		dstype xdg1 = xdg[0*ng+i];
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype nlg1 = nlg[0*ng+i];
		dstype nlg2 = nlg[1*ng+i];
		f[0*ng+i] = nlg1*param1*udg1*xdg2 - nlg2*param2*udg1*xdg1 + tau1*udg1;
	});
}

void KokkosFbou(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	if (ib == 1)
		KokkosFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

