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
    dstype v0 = v[0*N+i];
    dstype mu0 = mu[0];
    dstype mu1 = mu[1];
    dstype mu2 = mu[2];
    dstype mu3 = mu[3];
    dstype mu9 = mu[9];

    dstype x0 = pow(uq1, 2);
    dstype x1 = pow(uq0, -1);
    dstype x2 = x1*uq2;
    dstype x3 = uq10 - x2*uq8;
    dstype x4 = x1*x3;
    dstype x5 = x1*uq1;
    dstype x6 = uq5 - x5*uq4;
    dstype x7 = x1*x6;
    dstype x8 = -1.0 + mu0;
    dstype x9 = pow(uq2, 2);
    dstype x10 = pow(uq0, -2);
    dstype x11 = 0.5*(x0*x10 + x9*x10);
    dstype x12 = uq3 - uq0*x11;
    dstype x13 = x8*x12;
    dstype x14 = x1*x13;
    dstype x15 = sqrt(pow(x8, 3)*pow(x12, 3)*pow(mu3, 6)*pow(mu0, 3)/pow(uq0, 3))*(110.4 + mu9)/(mu1*(110.4 + 1.0*x14*pow(mu3, 2)*mu9*mu0));
    dstype x16 = 0.666666666666667*x15;
    dstype x17 = (-x4 + 2*x7)*x16;
    dstype x18 = uq6 - x2*uq4;
    dstype x19 = uq9 - x5*uq8;
    dstype x20 = 1.0*x15;
    dstype x21 = (x1*x18 + x1*x19)*x20;
    dstype x22 = x21 + x5*uq2;
    dstype x23 = uq2*x10;
    dstype x24 = uq1*x10;
    dstype x25 = x8*uq0;
    dstype x26 = x20*x10*mu0/(x8*mu2);
    dstype x27 = x14 + x1*uq3;
    dstype x28 = (2*x4 - x7)*x16;

    f[0 * N + i] = uq1 + v0*uq4;
    f[1 * N + i] = x13 + x17 + v0*uq5 + x0*x1;
    f[2 * N + i] = x22 + v0*uq6;
    f[3 * N + i] = uq1*x27 + v0*uq7 + x2*x21 + x26*(-uq4*x13 + x25*(uq7 - uq0*(x23*x18 + x6*x24) - uq4*x11)) + x5*x17;
    f[4 * N + i] = uq2 + v0*uq8;
    f[5 * N + i] = x22 + v0*uq9;
    f[6 * N + i] = x13 + x28 + v0*uq10 + x1*x9;
    f[7 * N + i] = uq2*x27 + v0*uq11 + x2*x28 + x26*(-uq8*x13 + x25*(uq11 - uq0*(x24*x19 + x3*x23) - uq8*x11)) + x5*x21;
  });
}

