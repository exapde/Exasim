void HdgFext1(dstype* f, dstype* J1, dstype* J2, dstype* J3, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* uhat, const dstype* n, const dstype* uext, const dstype* tau, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szuhat, const int szn, const int szuext, const int sztau, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Fext", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uhat0 = uhat[0*N+i];
    dstype uext0 = uext[0*N+i];

    f[0 * N + i] = uext0 - uhat0;
    J1[0 * N + i] = 0;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J3[0 * N + i] = -1;
  });
}

void HdgFext(dstype* f, dstype* f_udg, dstype* f_wdg, dstype* f_uhg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg,
           const dstype* nlg, const dstype* uext, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time,
           const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd,
           const int ncx, const int nco, const int ncw) {
    if (ib == 1 )
        HdgFext1(f, f_udg, f_wdg, f_uhg, xdg, udg, odg, wdg, uhg, nlg, uext, tau, uinf, param, time, modelnumber,
                        ng, nc, ncu, nd, ncx, nco, ncw, nc, ncu, nd, ncx);
}
