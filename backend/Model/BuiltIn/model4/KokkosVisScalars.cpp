void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("VisScalars", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype param1 = param[0];
		dstype param4 = param[3];
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype udg5 = udg[4*ng+i];
		dstype wdg1 = wdg[0*ng+i];
		dstype t2 = param1*udg1;
		dstype t3 = param1*udg2;
		dstype t4 = param1*udg3;
		dstype t5 = param1*udg4;
		dstype t6 = param1*udg5;
		f[0*ng+i] = t2;
		f[1*ng+i] = t3;
		f[2*ng+i] = t4;
		f[3*ng+i] = t5;
		f[4*ng+i] = t6;
		f[5*ng+i] = t2+t3+t4+t5+t6;
		f[6*ng+i] = param4*wdg1;
		f[7*ng+i] = param4*wdg1*(t2*7.139440410660612E+1+t3*6.250234383789392E+1+t4*3.332655693342354E+1+t5*3.569720205330306E+1+t6*3.125117191894696E+1)*8.314471468617452;
	});
}

