void KokkosMonitor(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)
{
	Kokkos::parallel_for("Monitor", ng, KOKKOS_LAMBDA(const size_t i) {
		int j = i%npe;
		int k = i/npe;
		dstype udg1 = udg[j+npe*0+npe*nc*k];
		f[j+npe*0+npe*nce*k] = udg1;
	});
}

