void KokkosInitu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)
{

  Kokkos::parallel_for("Initu", N, KOKKOS_LAMBDA(const size_t i) {
    int j = i%npe; 
    int k = i/npe; 


    f[j+npe*0 +npe*nce*k] = 0.0;
  });
}

