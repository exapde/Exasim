void HdgFbouonly1(dstype* f, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("FbouHdg", N, KOKKOS_LAMBDA(const size_t i) {
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

    dstype x0 = -1.0 + mu0;
    dstype x1 = pow(uq0, -2);
    dstype x2 = 0.5*(x1*pow(uq1, 2) + x1*pow(uq2, 2));
    dstype x3 = x0*(uq3 - x2*uq0);
    dstype x4 = pow(uq0, -1);
    dstype x5 = x4*uq2;
    dstype x6 = uq6 - x5*uq4;
    dstype x7 = x1*uq2;
    dstype x8 = x4*uq1;
    dstype x9 = uq5 - x8*uq4;
    dstype x10 = x1*uq1;
    dstype x11 = x0*uq0;
    dstype x12 = pow(mu1, -1);
    dstype x13 = x1*x12*mu0/(x0*mu2);
    dstype x14 = uq10 - x5*uq8;
    dstype x15 = x4*x14;
    dstype x16 = x4*x9;
    dstype x17 = 0.666666666666667*x12;
    dstype x18 = uq9 - x8*uq8;
    dstype x19 = x12*(x4*x18 + x4*x6);
    dstype x20 = x4*uq3 + x4*x3;

    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1;
    f[2 * N + i] = -uhat2;
    f[3 * N + i] = n0*(uq1*x20 + x13*(x11*(uq7 - uq0*(x6*x7 + x9*x10) - x2*uq4) - x3*uq4) + x5*x19 + (-x15 + 2*x16)*x8*x17) + n1*(uq2*x20 + x13*(x11*(uq11 - uq0*(x10*x18 + x7*x14) - x2*uq8) - x3*uq8) + x8*x19 + (2*x15 - x16)*x5*x17) + tau0*(-uhat3 + uq3);
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
    dstype n0 = n[0*N+i];
    dstype n1 = n[1*N+i];
    dstype mu0 = mu[0];
    dstype mu4 = mu[4];
    dstype mu5 = mu[5];
    dstype mu6 = mu[6];
    dstype mu7 = mu[7];

    dstype x0 = -mu5 + uq1;
    dstype x1 = pow(uhat0, -1);
    dstype x2 = pow(uhat0, -2);
    dstype x3 = x1*(uhat3 - 0.5*(x2*pow(uhat1, 2) + x2*pow(uhat2, 2))*uhat0)*(-1.0 + mu0);
    dstype x4 = sqrt(x3*mu0);
    dstype x5 = n0*uhat1 + n1*uhat2;
    dstype x6 = x1*x5;
    dstype x7 = x4 + x6;
    dstype x8 = tanh(100*x7);
    dstype x9 = x4*x6;
    dstype x10 = x3 + x1*uhat3;
    dstype x11 = x10 - x9;
    dstype x12 = -x11;
    dstype x13 = -x4 + x6;
    dstype x14 = -x13;
    dstype x15 = x14 + x6;
    dstype x16 = x15 + x4;
    dstype x17 = pow(x15, -1);
    dstype x18 = n0*uhat2 - n1*uhat1;
    dstype x19 = x2*pow(x18, 2);
    dstype x20 = 0.5*(x19 + x2*pow(x5, 2));
    dstype x21 = (x12 + x20)*x17;
    dstype x22 = x21*x16;
    dstype x23 = x10 + x9;
    dstype x24 = pow(x12 - x22 + x23, -1);
    dstype x25 = x1*x18;
    dstype x26 = x24*x25;
    dstype x27 = n1*x26;
    dstype x28 = x24*x21;
    dstype x29 = x8*(x27 - n0*x28);
    dstype x30 = tanh(100*x6);
    dstype x31 = x17*(1 + x24*x22);
    dstype x32 = x17*x16;
    dstype x33 = x30*(n0*x31 - x32*x27);
    dstype x34 = tanh(100*x13);
    dstype x35 = x26 - x32*x26;
    dstype x36 = x28 - x31;
    dstype x37 = (n0*x36 - n1*x35)*x34;
    dstype x38 = -mu6 + uq2;
    dstype x39 = n0*x26;
    dstype x40 = x8*(-x39 - n1*x28);
    dstype x41 = x30*(n1*x31 + x32*x39);
    dstype x42 = (n0*x35 + n1*x36)*x34;
    dstype x43 = -mu4 + uq0;
    dstype x44 = x24*(x12 + x19 + x21*x13);
    dstype x45 = x17*(x14 - x44*x16);
    dstype x46 = x45*x30;
    dstype x47 = x8*x44;
    dstype x48 = (1 - x44 - x45)*x34;
    dstype x49 = -mu7 + uq3;
    dstype x50 = x8*x24;
    dstype x51 = x32*x24;
    dstype x52 = x51*x30;
    dstype x53 = (-x24 + x51)*x34;
    dstype x54 = pow(n0, -1);
    dstype x55 = pow(n1, 2)*x54;
    dstype x56 = pow(n0 + x55, -1);
    dstype x57 = x56*x55;
    dstype x58 = x54*(1 - x57);
    dstype x59 = x56*x25;
    dstype x60 = n1*x54;
    dstype x61 = -x60*x59;
    dstype x62 = x61 + x7*x58;
    dstype x63 = x61 + x58*x13;
    dstype x64 = x61 + x6*x58;
    dstype x65 = x59*x30;
    dstype x66 = x56*x30;
    dstype x67 = -n1*x66;
    dstype x68 = x60*x56;
    dstype x69 = x59 + x7*x68;
    dstype x70 = x59 + x6*x68;
    dstype x71 = x59 + x68*x13;
    dstype x72 = x30*x25;

    f[0 * N + i] = -uhat0 + 0.5*(mu4 + uq0 + x0*(x29 + x33 + x37) + x38*(x40 + x41 + x42) + x43*(x46 + x47 + x48) + x49*(x50 - x52 + x53));
    f[1 * N + i] = -uhat1 + 0.5*(mu5 + uq1 + x0*(x57*x30 + x62*x29 + x63*x37 + x64*x33) + x38*(x67 + x62*x40 + x63*x42 + x64*x41) + x49*(x62*x50 + x63*x53 - x64*x52) + (x60*x65 + x62*x47 + x63*x48 + x64*x46)*x43);
    f[2 * N + i] = -uhat2 + 0.5*(mu6 + uq2 + x0*(x67 + x69*x29 + x70*x33 + x71*x37) + x38*(n0*x66 + x69*x40 + x70*x41 + x71*x42) + x43*(-x65 + x69*x47 + x70*x46 + x71*x48) + x49*(x69*x50 - x70*x52 + x71*x53));
    f[3 * N + i] = -uhat3 + 0.5*(mu7 + uq3 + x0*(-n1*x72 + x23*x29 + x33*x20 + x37*x11) + x38*(n0*x72 + x40*x23 + x41*x20 + x42*x11) + x49*(x50*x23 - x52*x20 + x53*x11) + (-x30*x19 + x46*x20 + x47*x23 + x48*x11)*x43);
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
}
