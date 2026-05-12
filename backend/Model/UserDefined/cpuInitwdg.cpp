void cpuInitwdg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int ncw_runtime, const int npe, const int ne)
{
    constexpr int nd = PdeModel::nd;
    constexpr int ncw = PdeModel::ncw;

    (void)modelnumber;
    (void)ng;

    if constexpr (ncw > 0) {
        for (int elem = 0; elem < ne; ++elem) {
            for (int j = 0; j < npe; ++j) {
                dstype x[nd];
                dstype w_local[ncw];

                for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
                PdeModel::cpuinitwdg(w_local, x, uinf, param);
                for (int k = 0; k < ncw; ++k) f[j + npe * k + npe * ncw_runtime * elem] = w_local[k];
            }
        }
    } else {
        (void)f; (void)xdg; (void)uinf; (void)param; (void)ncw_runtime; (void)npe; (void)ne;
    }
}
