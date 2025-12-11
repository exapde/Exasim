void KokkosSource(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Source", N, KOKKOS_LAMBDA(const size_t i) {
    dstype x1 = x[1*N+i];
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];
    dstype uq4 = uq[4*N+i];
    dstype uq5 = uq[5*N+i];
    dstype uq6 = uq[6*N+i];
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

    dstype x0 = pow(x1, -1);
    dstype x2 = pow(uq0, -1);
    dstype x3 = x2*uq2;
    dstype x4 = x2*uq1;
    dstype x5 = uq9 - x4*uq8;
    dstype x6 = -1.0 + mu0;
    dstype x7 = pow(uq2, 2);
    dstype x8 = pow(uq0, -2);
    dstype x9 = 0.5*(x8*pow(uq1, 2) + x8*x7);
    dstype x10 = uq3 - x9*uq0;
    dstype x11 = x6*x10;
    dstype x12 = x2*x11;
    dstype x13 = sqrt(pow(x6, 3)*pow(x10, 3)*pow(mu3, 6)*pow(mu0, 3)/pow(uq0, 3))*(110.4 + mu9)/(mu1*(110.4 + 1.0*x12*pow(mu3, 2)*mu9*mu0));
    dstype x14 = 1.0*x13;
    dstype x15 = uq10 - x3*uq8;
    dstype x16 = x2*x15;
    dstype x17 = x0*x3;
    dstype x18 = -x2*(uq5 - x4*uq4);
    dstype x19 = 0.666666666666667*x13;
    dstype x20 = x19*(2*x16 + x17 + x18);

    f[0 * N + i] = -x0*(uq2 + v0*uq8);
    f[1 * N + i] = -x0*(v0*uq9 + x14*(x2*x5 + x2*(uq6 - x3*uq4)) + x4*uq2);
    f[2 * N + i] = -x0*(x20 + v0*uq10 - x19*(-x16 - 2*x17 + x18) + x2*x7);
    f[3 * N + i] = -x0*(uq2*(x12 + x2*uq3) + v0*uq11 + x3*x20 + x8*x14*mu0*(-uq8*x11 + x6*uq0*(uq11 - uq0*(x5*x8*uq1 + x8*uq2*x15) - x9*uq8))/(x6*mu2));
  });
}

