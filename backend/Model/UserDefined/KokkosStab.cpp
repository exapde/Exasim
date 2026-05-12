template <class Model>
static void KokkosStabTemplate(dstype* f, const dstype* xdg,
                               const dstype* udg1, const dstype* udg2,
                               const dstype* odg1, const dstype* odg2,
                               const dstype* wdg1, const dstype* wdg2,
                               const dstype* uhg, const dstype* nlg,
                               const dstype* tau, const dstype* uinf,
                               const dstype* param, const dstype time,
                               const int modelnumber, const int ng,
                               const int nc_runtime,
                               const int ncu_runtime,
                               const int nd_runtime, const int ncx,
                               const int nco_runtime,
                               const int ncw_runtime)
{
    constexpr int nd = Model::nd;
    constexpr int ncu = Model::ncu;
    constexpr int nc = ncu * (1 + nd);
    constexpr int nco = Model::nco;
    constexpr int ncw = Model::ncw;

    (void)modelnumber;
    (void)nc_runtime;
    (void)ncu_runtime;
    (void)nd_runtime;
    (void)ncx;
    (void)nco_runtime;
    (void)ncw_runtime;

    Kokkos::parallel_for("Stab", ng, KOKKOS_LAMBDA(const size_t i) {
        dstype x[nd];
        dstype uq1[nc];
        dstype uq2[nc];
        dstype v1[(nco > 0) ? nco : 1];
        dstype v2[(nco > 0) ? nco : 1];
        dstype w1[(ncw > 0) ? ncw : 1];
        dstype w2[(ncw > 0) ? ncw : 1];
        dstype uh[ncu];
        dstype n[nd];
        dstype tau_local[ncu];
        dstype stab_local[ncu];

        for (int k = 0; k < nd; ++k) x[k] = xdg[k * ng + i];
        for (int k = 0; k < nc; ++k) {
            uq1[k] = udg1[k * ng + i];
            uq2[k] = udg2[k * ng + i];
        }
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) {
                v1[k] = odg1[k * ng + i];
                v2[k] = odg2[k * ng + i];
            }
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) {
                w1[k] = wdg1[k * ng + i];
                w2[k] = wdg2[k * ng + i];
            }
        }
        for (int k = 0; k < ncu; ++k) {
            uh[k] = uhg[k * ng + i];
            tau_local[k] = tau[k];
        }
        for (int k = 0; k < nd; ++k) n[k] = nlg[k * ng + i];

        Model::stab(stab_local, x, uq1, uq2, v1, v2, w1, w2, uh, n, tau_local,
                    param, uinf, time);

        for (int k = 0; k < ncu; ++k) f[k * ng + i] = stab_local[k];
    });
}

void KokkosStab(dstype* f, const dstype* xdg, const dstype* udg1, const dstype* udg2,  const dstype* odg1, const dstype* odg2,  const dstype* wdg1, const dstype* wdg2,  const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
    KokkosStabTemplate<PdeModel>(f, xdg, udg1, udg2, odg1, odg2, wdg1, wdg2,
                                 uhg, nlg, tau, uinf, param, time, modelnumber,
                                 ng, nc, ncu, nd, ncx, nco, ncw);
}
