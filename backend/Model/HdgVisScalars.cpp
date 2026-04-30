void HdgVisScalars(dstype* f, dstype* J1, dstype* J2, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("VisScalars", N, KOKKOS_LAMBDA(const size_t i) {
    dstype uq0 = uq[0*N+i];
    dstype uq1 = uq[1*N+i];
    dstype uq2 = uq[2*N+i];

    f[0 * N + i] = uq0;
    f[1 * N + i] = uq1 + uq2;
    J1[0 * N + i] = 1;
    J1[1 * N + i] = 0;
    J1[2 * N + i] = 0;
    J1[3 * N + i] = 1;
    J1[4 * N + i] = 0;
    J1[5 * N + i] = 1;
  });
}

