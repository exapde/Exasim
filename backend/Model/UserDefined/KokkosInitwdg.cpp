template <class Model>
static void KokkosInitwdgTemplate(dstype* f, const dstype* xdg,
                                  const dstype* uinf, const dstype* param,
                                  const int modelnumber, const int ng,
                                  const int ncx, const int ncw_runtime,
                                  const int npe, const int ne)
{
    constexpr int nd = Model::nd;
    constexpr int ncw = Model::ncw;

    (void)modelnumber;
    (void)ne;

    Kokkos::parallel_for("Initwdg", ng, KOKKOS_LAMBDA(const size_t i) {
        constexpr int nd = Model::nd;
        constexpr int ncw = Model::ncw;
        const int j = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);
        dstype x[nd];
        dstype w_local[(ncw > 0) ? ncw : 1];

        for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
        Model::initwdg(w_local, x, uinf, param);
        for (int k = 0; k < ncw; ++k) f[j + npe * k + npe * ncw_runtime * elem] = w_local[k];
    });
}

void KokkosInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
    KokkosInitwdgTemplate<PdeModel>(f, xdg, uinf, param, modelnumber, ng, ncx,
                                    nce, npe, ne);
}
