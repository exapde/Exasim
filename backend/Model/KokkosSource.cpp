void KokkosSource(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Source", N, KOKKOS_LAMBDA(const size_t i) {


    f[0 * N + i] = 0;
    f[1 * N + i] = 0;
    f[2 * N + i] = 0;
    f[3 * N + i] = 0;
  });
}

