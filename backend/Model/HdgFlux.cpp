void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Flux", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype odg4 = odg[3*ng+i];
		dstype odg5 = odg[4*ng+i];
		dstype odg6 = odg[5*ng+i];
		{
		dstype t2 = odg1*param2;
		dstype t3 = time-1.0;
		dstype t4 = t2-1.0;
		dstype t5 = odg1*t3;
		dstype t8 = -1.0/(t5-time);
		f[0*ng+i] = odg4*t8*udg1+odg6*param1*t4*udg2;
		f[1*ng+i] = odg5*t8*udg1+odg6*param1*t4*udg3;
		}
		{
		dstype t2 = odg1*param2;
		dstype t3 = time-1.0;
		dstype t4 = t2-1.0;
		dstype t5 = odg1*t3;
		dstype t7 = odg6*param1*t4;
		dstype t9 = -1.0/(t5-time);
		f_udg[0*ng+i] = odg4*t9;
		f_udg[1*ng+i] = odg5*t9;
		f_udg[2*ng+i] = t7;
		f_udg[3*ng+i] = 0.0;
		f_udg[4*ng+i] = 0.0;
		f_udg[5*ng+i] = t7;
		}
	});
}

