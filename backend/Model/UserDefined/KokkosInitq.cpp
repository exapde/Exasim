template <class Model>
static void KokkosInitqTemplate(dstype* f, const dstype* xdg,
                                const dstype* uinf, const dstype* param,
                                const int modelnumber, const int ng,
                                const int ncx, const int nc_runtime,
                                const int npe, const int ne)
{
    constexpr int nd = Model::nd;
    constexpr int nq = Model::ncu * Model::nd;

    (void)modelnumber;
    (void)ne;

    Kokkos::parallel_for("Initq", ng, KOKKOS_LAMBDA(const size_t i) {
        const int j = static_cast<int>(i % npe);
        const int elem = static_cast<int>(i / npe);
        dstype x[nd];
        dstype q_local[nq];

        for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
        Model::initq(q_local, x, uinf, param);
        for (int k = 0; k < nq; ++k) f[j + npe * k + npe * nc_runtime * elem] = q_local[k];
    });
}

void KokkosInitq(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
    KokkosInitqTemplate<PdeModel>(f, xdg, uinf, param, modelnumber, ng, ncx,
                                  nce, npe, ne);
}
