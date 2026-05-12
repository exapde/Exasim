template <class Model>
static void KokkosInituTemplate(dstype* f, const dstype* xdg,
                                const dstype* uinf, const dstype* param,
                                const int modelnumber, const int ng,
                                const int ncx, const int nc_runtime,
                                const int npe, const int ne)
{
    constexpr int nd = Model::nd;
    constexpr int ncu = Model::ncu;

    (void)modelnumber;
    (void)ne;

    Kokkos::parallel_for("Initu", ng, KOKKOS_LAMBDA(const size_t i) {
        const int j = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);
        dstype x[nd];
        dstype u_local[ncu];

        for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
        Model::initu(u_local, x, uinf, param);
        for (int k = 0; k < ncu; ++k) f[j + npe * k + npe * nc_runtime * elem] = u_local[k];
    });
}

void KokkosInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
    KokkosInituTemplate<PdeModel>(f, xdg, uinf, param, modelnumber, ng, ncx,
                                  nce, npe, ne);
}
