void KokkosFlux(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param3 = param[2];
		dstype param5 = param[4];
		dstype param10 = param[9];
		dstype param11 = param[10];
		dstype param12 = param[11];
		dstype param14 = param[13];
		dstype param16 = param[15];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = 1.0/3.141592653589793;
		dstype t3 = udg2*1.0E+3;
		dstype t4 = atan(t3);
		dstype t5 = t2*t4;
		dstype t6 = t5+1.0/2.0;
		dstype t7 = t6*udg2;
		f[0*ng+i] = udg3*(param16*(t7+3.183097800805168E-4)-param14*(t7-9.996816902199195E-1));
		f[1*ng+i] = 1.0/(param1*param1)*param10*param12*udg4*exp(-param11/(param3*param5*udg1));
	});
}

