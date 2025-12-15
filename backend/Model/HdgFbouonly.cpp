void HdgFbouonly1(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uhat0 = uhat[0*N+i];
    dstype uhat1 = uhat[1*N+i];
    dstype uhat2 = uhat[2*N+i];
    dstype uhat3 = uhat[3*N+i];
    dstype mu4 = mu[4];
    dstype mu5 = mu[5];
    dstype mu6 = mu[6];
    dstype mu7 = mu[7];


    f[0 * N + i] = mu4 - uhat0;
    f[1 * N + i] = mu5 - uhat1;
    f[2 * N + i] = mu6 - uhat2;
    f[3 * N + i] = mu7 - uhat3;
  });
}

void HdgFbouonly2(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype uhat1 = uhat[1*N+i];
    dstype uhat2 = uhat[2*N+i];
    dstype uhat3 = uhat[3*N+i];


    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1 + uq1;
    f[2 * N + i] = -uhat2 + uq2;
    f[3 * N + i] = -uhat3 + uq3;
  });
}

void HdgFbouonly3(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype uhat1 = uhat[1*N+i];
    dstype uhat2 = uhat[2*N+i];
    dstype uhat3 = uhat[3*N+i];
    dstype mu8 = mu[8];
    dstype mu9 = mu[9];
    dstype mu10 = mu[10];


    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1;
    f[2 * N + i] = -uhat2;
    f[3 * N + i] = -uhat3 + mu8*mu10*uhat0/mu9;
  });
}

void HdgFbouonly4(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
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
    dstype uhat0 = uhat[0*N+i];
    dstype uhat1 = uhat[1*N+i];
    dstype uhat2 = uhat[2*N+i];
    dstype uhat3 = uhat[3*N+i];
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype tau0 = tau[0];
    dstype mu0 = mu[0];
    dstype mu1 = mu[1];
    dstype mu2 = mu[2];
    dstype mu3 = mu[3];
    dstype mu9 = mu[9];

    dstype x0 = pow(uhat0, -1);
    dstype x1 = x0*uhat1;
    dstype x2 = uq5 - x1*uq4;
    dstype x3 = x0*x2;
    dstype x4 = x0*uhat2;
    dstype x5 = uq10 - x4*uq8;
    dstype x6 = x0*x5;
    dstype x7 = pow(uhat0, -2);
    dstype x8 = 0.5*(x7*pow(uhat1, 2) + x7*pow(uhat2, 2));
    dstype x9 = uhat3 - x8*uhat0;
    dstype x10 = -1.0 + mu0;
    dstype x11 = x9*x10;
    dstype x12 = x0*x11;
    dstype x13 = 1.0*mu0;
    dstype x14 = sqrt(pow(x9, 3)*pow(x10, 3)*pow(mu3, 6)*pow(mu0, 3)/pow(uhat0, 3))*(110.4 + mu9)/(mu1*(110.4 + x13*x12*pow(mu3, 2)*mu9));
    dstype x15 = 0.666666666666667*x14;
    dstype x16 = uq6 - x4*uq4;
    dstype x17 = uq9 - x1*uq8;
    dstype x18 = 1.0*(x0*x16 + x0*x17)*x14;
    dstype x19 = x12 + x0*uhat3;
    dstype x20 = x7*uhat2;
    dstype x21 = x7*uhat1;
    dstype x22 = x10*uhat0;
    dstype x23 = x7*x14*x13/(x10*mu2);

    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1;
    f[2 * N + i] = -uhat2;
    f[3 * N + i] = n0*(v0*uq7 + x19*uhat1 + x23*(-uq4*x11 + x22*(uq7 - uhat0*(x2*x21 + x20*x16) - x8*uq4)) + x4*x18 + (2*x3 - x6)*x1*x15) + n1*(v0*uq11 + x1*x18 + x19*uhat2 + x23*(-uq8*x11 + x22*(uq11 - uhat0*(x21*x17 + x5*x20) - x8*uq8)) + (-x3 + 2*x6)*x4*x15) + tau0*(-uhat3 + uq3);
  });
}

void HdgFbouonly5(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype v1 = v[1*N+i];
    dstype uhat0 = uhat[0*N+i];
    dstype uhat1 = uhat[1*N+i];
    dstype uhat2 = uhat[2*N+i];
    dstype uhat3 = uhat[3*N+i];


    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1;
    f[2 * N + i] = -uhat2;
    f[3 * N + i] = -uhat3 + v1*uhat0;
  });
}

void HdgFbouonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,
           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,
           const int ncx, const int nco, const int ncw) {
    if (ib == 1 )
        HdgFbouonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 2 )
        HdgFbouonly2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 3 )
        HdgFbouonly3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 4 )
        HdgFbouonly4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 5 )
        HdgFbouonly5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
}
