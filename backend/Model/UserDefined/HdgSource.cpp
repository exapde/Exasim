template <class Model>
static void HdgSourceTemplate(dstype* f, dstype* f_udg, dstype* f_wdg,
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

    Kokkos::parallel_for("HdgSource", ng, KOKKOS_LAMBDA(const size_t i) {
        constexpr int nd = Model::nd;
        constexpr int ncu = Model::ncu;
        constexpr int nc = ncu * (1 + nd);
        constexpr int nco = Model::nco;
        constexpr int ncw = Model::ncw;
        dstype x[nd];
        dstype uq[nc];
        dstype v[(nco > 0) ? nco : 1];
        dstype w[(ncw > 0) ? ncw : 1];
        dstype s_local[ncu];
        dstype s_uq[ncu * nc];

        for (int k = 0; k < nd; ++k) x[k] = xdg[k * ng + i];
        for (int k = 0; k < nc; ++k) uq[k] = udg[k * ng + i];
        for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        Model::source(s_local, x, uq, v, w, param, uinf, time);
        for (int k = 0; k < ncu; ++k) f[k * ng + i] = s_local[k];

        Model::source_jac_uq(s_uq, x, uq, v, w, param, uinf, time);
        for (int k = 0; k < ncu * nc; ++k) f_udg[k * ng + i] = s_uq[k];

        dstype s_w[(ncu * ncw > 0) ? ncu * ncw : 1];
        Model::source_jac_w(s_w, x, uq, v, w, param, uinf, time);
        for (int k = 0; k < ncu * ncw; ++k) f_wdg[k * ng + i] = s_w[k];
    });
}

void HdgSource(dstype* f, dstype* f_udg, dstype* f_wdg, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw)
{
    HdgSourceTemplate<PdeModel>(f, f_udg, f_wdg, xdg, udg, odg, wdg, uinf,
                                param, time, modelnumber, ng, nc, ncu, nd,
                                ncx, nco, ncw);
}
