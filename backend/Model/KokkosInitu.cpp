void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	Kokkos::parallel_for("Initu", ng, KOKKOS_LAMBDA(const size_t i) {
		int j = i%npe;
		int k = i/npe;
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		f[j+npe*0+npe*nce*k] = xdg1;
	});
}

