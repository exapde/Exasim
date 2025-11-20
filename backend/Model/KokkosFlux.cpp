void KokkosFlux(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Flux", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];
    dstype uq4 = uq[4*N+i];
    dstype uq5 = uq[5*N+i];
    dstype uq6 = uq[6*N+i];
    dstype uq7 = uq[7*N+i];
    dstype uq8 = uq[8*N+i];
    dstype uq9 = uq[9*N+i];
    dstype uq10 = uq[10*N+i];
    dstype uq11 = uq[11*N+i];
    dstype mu0 = mu[0];
    dstype mu1 = mu[1];
    dstype mu2 = mu[2];

    dstype x0 = -1.0 + mu0;
    dstype x1 = pow(uq2, 2);
    dstype x2 = pow(uq0, -2);
    dstype x3 = pow(uq1, 2);
    dstype x4 = 0.5*(x2*x1 + x2*x3);
    dstype x5 = x0*(uq3 - x4*uq0);
    dstype x6 = pow(uq0, -1);
    dstype x7 = x6*uq2;
    dstype x8 = uq10 - x7*uq8;
    dstype x9 = x6*x8;
    dstype x10 = x6*uq1;
    dstype x11 = uq5 - uq4*x10;
    dstype x12 = x6*x11;
    dstype x13 = pow(mu1, -1);
    dstype x14 = 0.666666666666667*x13;
    dstype x15 = x14*(2*x12 - x9);
    dstype x16 = uq6 - x7*uq4;
    dstype x17 = uq9 - uq8*x10;
    dstype x18 = (x6*x16 + x6*x17)*x13;
    dstype x19 = x18 + x7*uq1;
    dstype x20 = x2*uq2;
    dstype x21 = x2*uq1;
    dstype x22 = x0*uq0;
    dstype x23 = x2*x13*mu0/(x0*mu2);
    dstype x24 = x6*uq3 + x6*x5;
    dstype x25 = x14*(-x12 + 2*x9);

    f[0 * N + i] = uq1;
    f[1 * N + i] = x15 + x5 + x3*x6;
    f[2 * N + i] = x19;
    f[3 * N + i] = uq1*x24 + x15*x10 + x23*(x22*(uq7 - x4*uq4 - (x20*x16 + x21*x11)*uq0) - x5*uq4) + x7*x18;
    f[4 * N + i] = uq2;
    f[5 * N + i] = x19;
    f[6 * N + i] = x25 + x5 + x1*x6;
    f[7 * N + i] = uq2*x24 + x10*x18 + x23*(x22*(uq11 - uq0*(x21*x17 + x8*x20) - x4*uq8) - x5*uq8) + x7*x25;
  });
}

