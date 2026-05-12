template <class Model>
static void HdgFluxTemplate(dstype* f, dstype* f_udg, dstype* f_wdg,
                            const dstype* xdg, const dstype* udg,
                            const dstype* odg, const dstype* wdg,
                            const dstype* uinf, const dstype* param,
                            const dstype time, const int modelnumber,
                            const int ng, const int nc_runtime,
                            const int ncu_runtime, const int nd_runtime,
                            const int ncx, const int nco_runtime,
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

    Kokkos::parallel_for("HdgFlux", ng, KOKKOS_LAMBDA(const size_t i) {
        dstype x[nd];
        dstype uq[nc];
        dstype v[(nco > 0) ? nco : 1];
        dstype w[(ncw > 0) ? ncw : 1];
        dstype f_local[ncu * nd];
        dstype f_uq[ncu * nd * nc];

        for (int k = 0; k < nd; ++k) x[k] = xdg[k * ng + i];
        for (int k = 0; k < nc; ++k) uq[k] = udg[k * ng + i];
        if constexpr (nco > 0) {
            for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        }
        if constexpr (ncw > 0) {
            for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];
        }

        Model::flux(f_local, x, uq, v, w, param, uinf, time);
        for (int k = 0; k < ncu * nd; ++k) f[k * ng + i] = f_local[k];

        Model::flux_jac_uq(f_uq, x, uq, v, w, param, uinf, time);
        for (int k = 0; k < ncu * nd * nc; ++k) f_udg[k * ng + i] = f_uq[k];

        if constexpr (ncw > 0) {
            dstype f_w[ncu * nd * ncw];
            Model::flux_jac_w(f_w, x, uq, v, w, param, uinf, time);
            for (int k = 0; k < ncu * nd * ncw; ++k) f_wdg[k * ng + i] = f_w[k];
        }
    });
}

void HdgFlux(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
    HdgFluxTemplate<PdeModel>(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf,
                              param, time, modelnumber, ng, nc, ncu, nd, ncx,
                              nco, ncw);
}
