void KokkosInitu(dstype* f, const dstype* x, const dstype* eta, const dstype* mu, const int modelnumber, const int N, const int ncx, const int nce, const int npe, const int ne)
{

  Kokkos::parallel_for("Initu", N, KOKKOS_LAMBDA(const size_t i) {
    int p = i%npe; 
    int e = i/npe; 
    dstype mu4 = mu[4];
    dstype mu5 = mu[5];
    dstype mu6 = mu[6];
    dstype mu7 = mu[7];


    f[p+npe*0 +npe*nce*e] = mu4;
    f[p+npe*1 +npe*nce*e] = mu5;
    f[p+npe*2 +npe*nce*e] = mu6;
    f[p+npe*3 +npe*nce*e] = mu7;
  });
}

