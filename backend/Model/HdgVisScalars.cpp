void HdgVisScalars(dstype* f, dstype* J1, dstype* J2, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("VisScalars", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];
    dstype uq3 = uq[3*N+i];
    dstype x0 = pow(uq0, -1);
    dstype x1 = x0*uq1;
    dstype x2 = x0*uq2;
    dstype x3 = pow(uq2, 2);
    dstype x4 = pow(uq1, 2);
    dstype x5 = pow(uq0, -2);
    dstype x6 = 0.5*x5;

    f[0 * N + i] = uq0;
    f[1 * N + i] = x1;
    f[2 * N + i] = x2;
    f[3 * N + i] = 0.4*(uq3 - 0.5*(x0*x3 + x0*x4));
    J1[0 * N + i] = 1;
    J1[1 * N + i] = -x5*uq1;
    J1[2 * N + i] = -x5*uq2;
    J1[3 * N + i] = 0.4*(x3*x6 + x4*x6);
    J1[4 * N + i] = 0;
    J1[5 * N + i] = x0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = -0.4*x1;
    J1[8 * N + i] = 0;
    J1[9 * N + i] = 0;
    J1[10 * N + i] = x0;
    J1[11 * N + i] = -0.4*x2;
    J1[12 * N + i] = 0;
    J1[13 * N + i] = 0;
    J1[14 * N + i] = 0;
    J1[15 * N + i] = 0.4;
    J1[16 * N + i] = 0;
    J1[17 * N + i] = 0;
    J1[18 * N + i] = 0;
    J1[19 * N + i] = 0;
    J1[20 * N + i] = 0;
    J1[21 * N + i] = 0;
    J1[22 * N + i] = 0;
    J1[23 * N + i] = 0;
    J1[24 * N + i] = 0;
    J1[25 * N + i] = 0;
    J1[26 * N + i] = 0;
    J1[27 * N + i] = 0;
    J1[28 * N + i] = 0;
    J1[29 * N + i] = 0;
    J1[30 * N + i] = 0;
    J1[31 * N + i] = 0;
    J1[32 * N + i] = 0;
    J1[33 * N + i] = 0;
    J1[34 * N + i] = 0;
    J1[35 * N + i] = 0;
    J1[36 * N + i] = 0;
    J1[37 * N + i] = 0;
    J1[38 * N + i] = 0;
    J1[39 * N + i] = 0;
    J1[40 * N + i] = 0;
    J1[41 * N + i] = 0;
    J1[42 * N + i] = 0;
    J1[43 * N + i] = 0;
    J1[44 * N + i] = 0;
    J1[45 * N + i] = 0;
    J1[46 * N + i] = 0;
    J1[47 * N + i] = 0;
  });
}

