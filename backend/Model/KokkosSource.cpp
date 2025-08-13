void KokkosSource(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Source", N, KOKKOS_LAMBDA(const size_t i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];


    f[0 * N + i] = 2*pow(acos(-1), 2)*sin(acos(-1)*x1)*sin(acos(-1)*x0);
  });
}

