void KokkosVisScalars(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
	Kokkos::parallel_for("VisScalars", ng, KOKKOS_LAMBDA(const size_t i) {
		dstype udg1 = udg[0*ng+i];
		dstype udg2 = udg[1*ng+i];
		dstype udg3 = udg[2*ng+i];
		dstype udg4 = udg[3*ng+i];
		dstype t2 = 1.0/udg1;
		f[0*ng+i] = udg1;
		f[1*ng+i] = t2*udg2;
		f[2*ng+i] = t2*udg3;
		f[3*ng+i] = udg4*(2.0/5.0)-(t2*(udg2*udg2))/5.0-(t2*(udg3*udg3))/5.0;
	});
}

