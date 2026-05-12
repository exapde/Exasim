template <class Model>
static void KokkosQoIboundaryTemplate(dstype* f, const dstype* xdg,
                                      const dstype* udg, const dstype* odg,
                                      const dstype* wdg, const dstype* uhg,
                                      const dstype* nlg, const dstype* tau,
                                      const dstype* uinf, const dstype* param,
                                      const dstype time, const int modelnumber,
                                      const int ib, const int ng,
                                      const int nvqoi_runtime,
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
    constexpr int kMax = 16;

    (void)modelnumber;
    (void)ncu_runtime;
    (void)nd_runtime;
    (void)ncx;
    (void)nco_runtime;
    (void)ncw_runtime;

    Kokkos::parallel_for("QoIboundary", ng, KOKKOS_LAMBDA(const size_t i) {
        dstype x[nd];
        dstype uq[nc];
        dstype v[(nco > 0) ? nco : 1];
        dstype w[(ncw > 0) ? ncw : 1];
        dstype uh[ncu];
        dstype n[nd];
        dstype tau_local[ncu];
        dstype out_local[kMax];

        for (int k = 0; k < kMax; ++k) out_local[k] = 0.0;
        for (int k = 0; k < nd; ++k) x[k] = xdg[k * ng + i];
        for (int k = 0; k < nc; ++k) uq[k] = udg[k * ng + i];
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        }
        for (int k = 0; k < ncu; ++k) {
            uh[k] = uhg[k * ng + i];
            tau_local[k] = tau[k];
        }
        for (int k = 0; k < nd; ++k) n[k] = nlg[k * ng + i];

        Model::qoi_boundary(out_local, ib, x, uq, v, w, uh, n, tau_local,
                            param, uinf, time);

        for (int k = 0; k < nvqoi_runtime; ++k) f[k * ng + i] = out_local[k];
    });
}

void KokkosQoIboundary(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uhg, const dstype* nlg, const dstype* tau, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ib, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
    KokkosQoIboundaryTemplate<PdeModel>(f, xdg, udg, odg, wdg, uhg, nlg, tau,
                                        uinf, param, time, modelnumber, ib, ng,
                                        nc, ncu, nd, ncx, nco, ncw);
}
