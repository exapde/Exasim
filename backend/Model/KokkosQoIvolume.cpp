void KokkosQoIvolume(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("QoIvolume", N, KOKKOS_LAMBDA(const size_t i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];
    dstype uq0 = uq[0*N+i];


    f[0 * N + i] = pow(uq0 - sin(acos(-1)*x1)*sin(acos(-1)*x0), 2);
    f[1 * N + i] = uq0;
  });
}

