void HdgSource(dstype* f, dstype* J1, dstype* J2, const dstype* x, const dstype* uq, const dstype* v, const dstype* w, const dstype* eta, const dstype* mu, const dstype t, const int modelnumber, const int N, const int szx, const int szuq, const int szv, const int szw, const int szeta, const int szmu)
{

  Kokkos::parallel_for("Source", N, KOKKOS_LAMBDA(const size_t i) {

    f[0 * N + i] = 0;
    f[1 * N + i] = 0;
    f[2 * N + i] = 0;
    f[3 * N + i] = 0;
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
  });
}

