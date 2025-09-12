void HdgQoIboundary(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("QoIboundary", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype tau0 = tau[0];
    dstype mu0 = mu[0];
    dstype x0 = n0*mu0;
    dstype x1 = n1*mu0;

    f[0 * N + i] = tau0*(-uhat0 + uq0) + x0*uq1 + x1*uq2;
    J1[0 * N + i] = tau0;
    J1[1 * N + i] = x0;
    J1[2 * N + i] = x1;
    J3[0 * N + i] = -tau0;
  });
}

