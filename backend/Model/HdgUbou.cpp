void HdgUbou1(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Ubou", N, KOKKOS_LAMBDA(const size_t i) {

    f[0 * N + i] = 0.0;
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J3[0 * N + i] = 0;
  });
}

