void HdgFbou1(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
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
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = 0;
    J1[8 * N + i] = 0;
    J1[9 * N + i] = 0;
    J1[10 * N + i] = 0;
    J1[11 * N + i] = 0;
    J1[12 * N + i] = 0;
    J1[13 * N + i] = 0;
    J1[14 * N + i] = 0;
    J1[15 * N + i] = 0;
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
    J3[0 * N + i] = -1;
    J3[1 * N + i] = 0;
    J3[2 * N + i] = 0;
    J3[3 * N + i] = 0;
    J3[4 * N + i] = 0;
    J3[5 * N + i] = -1;
    J3[6 * N + i] = 0;
    J3[7 * N + i] = 0;
    J3[8 * N + i] = 0;
    J3[9 * N + i] = 0;
    J3[10 * N + i] = -1;
    J3[11 * N + i] = 0;
    J3[12 * N + i] = 0;
    J3[13 * N + i] = 0;
    J3[14 * N + i] = 0;
    J3[15 * N + i] = -1;
  });
}

void HdgFbou2(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
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
    J1[0 * N + i] = 1;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 1;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = 0;
    J1[8 * N + i] = 0;
    J1[9 * N + i] = 0;
    J1[10 * N + i] = 1;
    J1[11 * N + i] = 0;
    J1[12 * N + i] = 0;
    J1[13 * N + i] = 0;
    J1[14 * N + i] = 0;
    J1[15 * N + i] = 1;
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
    J3[0 * N + i] = -1;
    J3[1 * N + i] = 0;
    J3[2 * N + i] = 0;
    J3[3 * N + i] = 0;
    J3[4 * N + i] = 0;
    J3[5 * N + i] = -1;
    J3[6 * N + i] = 0;
    J3[7 * N + i] = 0;
    J3[8 * N + i] = 0;
    J3[9 * N + i] = 0;
    J3[10 * N + i] = -1;
    J3[11 * N + i] = 0;
    J3[12 * N + i] = 0;
    J3[13 * N + i] = 0;
    J3[14 * N + i] = 0;
    J3[15 * N + i] = -1;
  });
}

void HdgFbou3(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int sztau, const int szeta, const int szmu)
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
    dstype x0 = mu8*mu10/mu9;

    f[0 * N + i] = -uhat0 + uq0;
    f[1 * N + i] = -uhat1;
    f[2 * N + i] = -uhat2;
    f[3 * N + i] = -uhat3 + x0*uhat0;
    J1[0 * N + i] = 1;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 0;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 0;
    J1[6 * N + i] = 0;
    J1[7 * N + i] = 0;
    J1[8 * N + i] = 0;
    J1[9 * N + i] = 0;
    J1[10 * N + i] = 0;
    J1[11 * N + i] = 0;
    J1[12 * N + i] = 0;
    J1[13 * N + i] = 0;
    J1[14 * N + i] = 0;
    J1[15 * N + i] = 0;
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
    J3[0 * N + i] = -1;
    J3[1 * N + i] = 0;
    J3[2 * N + i] = 0;
    J3[3 * N + i] = x0;
    J3[4 * N + i] = 0;
    J3[5 * N + i] = -1;
    J3[6 * N + i] = 0;
    J3[7 * N + i] = 0;
    J3[8 * N + i] = 0;
    J3[9 * N + i] = 0;
    J3[10 * N + i] = -1;
    J3[11 * N + i] = 0;
    J3[12 * N + i] = 0;
    J3[13 * N + i] = 0;
    J3[14 * N + i] = 0;
    J3[15 * N + i] = -1;
  });
}

void HdgFbou(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
           const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,
           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,
           const int ncx, const int nco, const int ncw) {
    if (ib == 1 )
        HdgFbou1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 2 )
        HdgFbou2(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
    else if (ib == 3 )
        HdgFbou3(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd);
}
