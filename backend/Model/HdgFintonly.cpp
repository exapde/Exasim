void HdgFintonly1(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Fint", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq3 = uq[3*N+i];
    dstype uq4 = uq[4*N+i];
    dstype uq5 = uq[5*N+i];
    dstype uq6 = uq[6*N+i];
    dstype uq7 = uq[7*N+i];
    dstype uq8 = uq[8*N+i];
    dstype uq9 = uq[9*N+i];
    dstype uq10 = uq[10*N+i];
    dstype uq11 = uq[11*N+i];
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

    dstype x0 = -1.0 + mu0;
    dstype x1 = pow(uhat0, -2);
    dstype x2 = 0.5*(x1*pow(uhat1, 2) + x1*pow(uhat2, 2));
    dstype x3 = uhat3 - x2*uhat0;
    dstype x4 = x0*x3;
    dstype x5 = pow(uhat0, -1);
    dstype x6 = x5*uq8;
    dstype x7 = x1*uhat2;
    dstype x8 = x1*uhat1;
    dstype x9 = x0*uhat0;
    dstype x10 = 1.0*mu0;
    dstype x11 = x1*x10*sqrt(pow(x0, 3)*pow(x3, 3)*pow(mu3, 6)*pow(mu0, 3)/pow(uhat0, 3))*(110.4 + mu9)/(x0*mu2*mu1*(110.4 + x4*x5*x10*pow(mu3, 2)*mu9));
    dstype x12 = x5*uq4;

    f[0 * N + i] = tau0*(-uhat3 + uq3) + n0*x11*(-x4*uq4 + x9*(uq7 - x2*uq4 - (x7*(uq6 - x12*uhat2) + x8*(uq5 - x12*uhat1))*uhat0)) + n1*x11*(-x4*uq8 + x9*(uq11 - uhat0*(x7*(uq10 - x6*uhat2) + x8*(uq9 - x6*uhat1)) - x2*uq8));
  });
}

void HdgFintonly(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,
           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,
           const int ncx, const int nco, const int ncw) {
    if (ib == 1 )
        HdgFintonly1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
}
