template <class Model>
static void KokkosInitodgTemplate(dstype* f, const dstype* xdg,
                                  const dstype* uinf, const dstype* param,
                                  const int modelnumber, const int ng,
                                  const int ncx, const int nco_runtime,
                                  const int npe, const int ne)
{
    constexpr int nd = Model::nd;
    constexpr int nco = Model::nco;

    (void)modelnumber;
    (void)ne;

    Kokkos::parallel_for("Initodg", ng, KOKKOS_LAMBDA(const size_t i) {
        const int j = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);

        for (int k = 0; k < nco_runtime; ++k) {
            f[j + npe * k + npe * nco_runtime * elem] = 0.0;
        }

        if constexpr (nco > 0) {
            dstype x[nd];
            dstype odg_local[nco];

            for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
            for (int k = 0; k < nco; ++k) odg_local[k] = 0.0;

            Model::initodg(odg_local, x, uinf, param);

            const int ncopy = (nco_runtime < nco) ? nco_runtime : nco;
            for (int k = 0; k < ncopy; ++k) {
                f[j + npe * k + npe * nco_runtime * elem] = odg_local[k];
            }
        }
    });
}

void KokkosInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
    KokkosInitodgTemplate<PdeModel>(f, xdg, uinf, param, modelnumber, ng, ncx,
                                    nce, npe, ne);
}
