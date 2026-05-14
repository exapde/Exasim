template <class Model>
static void KokkosEoSduTemplate(dstype* f, const dstype* xdg,
                                const dstype* udg, const dstype* odg,
                                const dstype* wdg, const dstype* uinf,
                                const dstype* param, const dstype time,
                                const int modelnumber, const int ng,
                                const int nc_runtime,
                                const int ncu_runtime,
                                const int nd_runtime, const int ncx,
                                const int nco_runtime,
                                const int ncw_runtime,
                                const int nce_runtime, const int npe,
                                const int ne)
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
    (void)nce_runtime;
    (void)npe;
    (void)ne;

    Kokkos::parallel_for("EoSdu", ng, KOKKOS_LAMBDA(const size_t i) {
        constexpr int nd = Model::nd;
        constexpr int ncu = Model::ncu;
        constexpr int nc = ncu * (1 + nd);
        constexpr int nco = Model::nco;
        constexpr int ncw = Model::ncw;
        dstype x[nd];
        dstype uq[nc];
        dstype v[(nco > 0) ? nco : 1];
        dstype w[(ncw > 0) ? ncw : 1];
        dstype eos_du_local[ncu * nc];

        for (int k = 0; k < nd; ++k) x[k] = xdg[k * ng + i];
        for (int k = 0; k < nc; ++k) uq[k] = udg[k * ng + i];
        for (int k = 0; k < nco; ++k) v[k] = odg[k * ng + i];
        for (int k = 0; k < ncw; ++k) w[k] = wdg[k * ng + i];

        Model::eos_du(eos_du_local, x, uq, v, w, param, uinf, time);

        for (int k = 0; k < ncu * nc; ++k) f[k * ng + i] = eos_du_local[k];
    });
}

void KokkosEoSdu(dstype* f, const dstype* xdg, const dstype* udg, const dstype* odg, const dstype* wdg, const dstype* uinf, const dstype* param, const dstype time, const int modelnumber, const int ng, const int nc, const int ncu, const int nd, const int ncx, const int nco, const int ncw, const int nce, const int npe, const int ne)
{
    KokkosEoSduTemplate<PdeModel>(f, xdg, udg, odg, wdg, uinf, param, time,
                                  modelnumber, ng, nc, ncu, nd, ncx, nco,
                                  ncw, nce, npe, ne);
}
