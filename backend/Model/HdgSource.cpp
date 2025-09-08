void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype odg1 = odg[0*ng+i];
		dstype odg2 = odg[1*ng+i];
		dstype odg3 = odg[2*ng+i];
		dstype odg4 = odg[3*ng+i];
		dstype odg5 = odg[4*ng+i];
		{
		dstype t2 = time-1.0;
		dstype t3 = odg1*t2;
		f[0*ng+i] = (udg1*(-odg1+t2*((odg2*odg4)/(t3-time)+(odg3*odg5)/(t3-time))+1.0))/(t3-time);
		}
		{
		dstype t2 = time-1.0;
		dstype t3 = odg1*t2;
		f_udg[0*ng+i] = (-odg1+t2*((odg2*odg4)/(t3-time)+(odg3*odg5)/(t3-time))+1.0)/(t3-time);
		f_udg[1*ng+i] = 0.0;
		f_udg[2*ng+i] = 0.0;
		}
	});
}

