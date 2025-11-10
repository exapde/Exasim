void KokkosSource(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("Source", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype xdg2 = xdg[1*ng+i];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = udg3*udg3;
		dstype t3 = 1.0/udg1;
		dstype t5 = 1.0/xdg2;
		dstype t4 = t3*t3;
		f[0*ng+i] = -t5*udg3;
		f[1*ng+i] = -t3*t5*udg2*udg3;
		f[2*ng+i] = -t2*t3*t5;
		f[3*ng+i] = t5*udg3*((t2*t4)/5.0-t3*udg4*(7.0/5.0)+(t4*(udg2*udg2))/5.0);
	});
}

