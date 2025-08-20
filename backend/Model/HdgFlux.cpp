void HdgFlux(dstype* f, dstype* J1, dstype* J2, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Flux", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];
    dstype mu0 = mu[0];

    f[0 * N + i] = uq1*mu0;
    f[1 * N + i] = uq2*mu0;
    f[2 * N + i] = uq3*mu0;
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = mu0;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = mu0;
    J1[8 * N + i] = 0;
    J1[9 * N + i] = 0;
    J1[10 * N + i] = 0;
    J1[11 * N + i] = mu0;
  });
}

