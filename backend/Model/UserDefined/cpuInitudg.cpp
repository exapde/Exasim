void cpuInitudg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nc_runtime, const int npe, const int ne)
{
    constexpr int nd = PdeModel::nd;
    constexpr int nqu = PdeModel::ncu * (1 + PdeModel::nd);

    (void)modelnumber;
    (void)ng;

    for (int elem = 0; elem < ne; ++elem) {
        for (int j = 0; j < npe; ++j) {
            dstype x[nd];
            dstype udg_local[nqu];

            for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
            PdeModel::cpuinitudg(udg_local, x, uinf, param);
            for (int k = 0; k < nqu; ++k) f[j + npe * k + npe * nc_runtime * elem] = udg_local[k];
        }
    }
}
