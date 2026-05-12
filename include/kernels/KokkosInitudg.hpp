template <class Model>
static void KokkosInitudgTemplate(dstype* f, const dstype* xdg,
                                  const dstype* uinf, const dstype* param,
                                  const int modelnumber, const int ng,
                                  const int ncx, const int nc_runtime,
                                  const int npe, const int ne)
{
    constexpr int nd = Model::nd;
    constexpr int nqu = Model::ncu * (1 + Model::nd);

    (void)modelnumber;
    (void)ne;

    Kokkos::parallel_for("Initudg", ng, KOKKOS_LAMBDA(const size_t i) {
        const int j = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);
        dstype x[nd];
        dstype udg_local[nqu];

        for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
        Model::initudg(udg_local, x, uinf, param);
        for (int k = 0; k < nqu; ++k) f[j + npe * k + npe * nc_runtime * elem] = udg_local[k];
    });
}

void KokkosInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
    KokkosInitudgTemplate<PdeModel>(f, xdg, uinf, param, modelnumber, ng, ncx,
                                    nce, npe, ne);
}
