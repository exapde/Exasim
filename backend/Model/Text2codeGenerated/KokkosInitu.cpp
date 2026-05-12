void KokkosInitu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)
{

  Kokkos::parallel_for("Initu", N, KOKKOS_LAMBDA(const size_t i) {
    int p = i%npe; 
    int e = i/npe; 


    f[p+npe*0 +npe*nce*e] = 0.0;
  });
}

