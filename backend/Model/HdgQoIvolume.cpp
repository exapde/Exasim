void HdgQoIvolume(dstype* f, dstype* J1, dstype* J2, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("QoIvolume", N, KOKKOS_LAMBDA(const size_t i) {
    dstype x0 = x[0*N+i];
    dstype x1 = x[1*N+i];
    dstype uq0 = uq[0*N+i];
    dstype x2 = uq0 - sin(acos(-1)*x1)*sin(acos(-1)*x0);

    f[0 * N + i] = pow(x2, 2);
    f[1 * N + i] = uq0;
    J1[0 * N + i] = 2*x2;
    J1[1 * N + i] = 1;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 0;
  });
}

