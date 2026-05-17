void KokkosVisScalars(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("VisScalars", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];

    dstype x0 = pow(uq0, -1);

    f[0 * N + i] = uq0;
    f[1 * N + i] = x0*uq1;
    f[2 * N + i] = x0*uq2;
    f[3 * N + i] = 0.4*(uq3 - 0.5*(x0*pow(uq1, 2) + x0*pow(uq2, 2)));
  });
}

