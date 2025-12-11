void KokkosFlux(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Flux", N, KOKKOS_LAMBDA(const size_t i) {
    dstype x1 = x[1*N+i];
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
    dstype v0 = v[0*N+i];
    dstype mu0 = mu[0];
    dstype mu1 = mu[1];
    dstype mu2 = mu[2];
    dstype mu3 = mu[3];
    dstype mu9 = mu[9];

    dstype x0 = pow(uq1, 2);
    dstype x2 = pow(uq0, -1);
    dstype x3 = -1.0 + mu0;
    dstype x4 = pow(uq2, 2);
    dstype x5 = pow(uq0, -2);
    dstype x6 = 0.5*(x0*x5 + x4*x5);
    dstype x7 = uq3 - x6*uq0;
    dstype x8 = x3*x7;
    dstype x9 = x2*uq2;
    dstype x10 = x9/x1;
    dstype x11 = x2*uq1;
    dstype x12 = uq5 - uq4*x11;
    dstype x13 = x2*x12;
    dstype x14 = uq10 - x9*uq8;
    dstype x15 = x2*x14;
    dstype x16 = x2*x8;
    dstype x17 = (110.4 + mu9)*sqrt(pow(x3, 3)*pow(x7, 3)*pow(mu3, 6)*pow(mu0, 3)/pow(uq0, 3))/(mu1*(110.4 + 1.0*x16*pow(mu3, 2)*mu9*mu0));
    dstype x18 = 0.666666666666667*x17;
    dstype x19 = x18*(x10 + 2*x13 - x15);
    dstype x20 = uq6 - x9*uq4;
    dstype x21 = uq9 - uq8*x11;
    dstype x22 = 1.0*x17;
    dstype x23 = (x2*x20 + x2*x21)*x22;
    dstype x24 = x23 + x9*uq1;
    dstype x25 = x16 + x2*uq3;
    dstype x26 = x5*uq2;
    dstype x27 = x5*uq1;
    dstype x28 = x3*uq0;
    dstype x29 = x5*x22*mu0/(x3*mu2);
    dstype x30 = x18*(x10 - x13 + 2*x15);

    f[0 * N + i] = uq1 + v0*uq4;
    f[1 * N + i] = x19 + x8 + v0*uq5 + x0*x2;
    f[2 * N + i] = x24 + v0*uq6;
    f[3 * N + i] = uq1*x25 + v0*uq7 + x11*x19 + x29*(x28*(uq7 - x6*uq4 - (x20*x26 + x27*x12)*uq0) - x8*uq4) + x9*x23;
    f[4 * N + i] = uq2 + v0*uq8;
    f[5 * N + i] = x24 + v0*uq9;
    f[6 * N + i] = x30 + x8 + v0*uq10 + x2*x4;
    f[7 * N + i] = uq2*x25 + v0*uq11 + x23*x11 + x29*(x28*(uq11 - x6*uq8 - (x21*x27 + x26*x14)*uq0) - x8*uq8) + x9*x30;
  });
}

