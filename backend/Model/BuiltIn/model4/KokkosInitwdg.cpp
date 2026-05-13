void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	Kokkos::parallel_for("Initwdg", ng, KOKKOS_LAMBDA(const size_t i) {
		int j = i%npe;
		int k = i/npe;
		f[j+npe*0+npe*nce*k] = 1.0;
	});
}

