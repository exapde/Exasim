void cpuInitodg(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nco_runtime, const int npe, const int ne)
{
    constexpr int nd = PdeModel::nd;
    constexpr int nco = PdeModel::nco;

    (void)modelnumber;
    (void)ng;

    for (int elem = 0; elem < ne; ++elem) {
        for (int j = 0; j < npe; ++j) {
            for (int k = 0; k < nco_runtime; ++k) {
                f[j + npe * k + npe * nco_runtime * elem] = 0.0;
            }

            dstype x[nd];
            dstype odg_local[(nco > 0) ? nco : 1];

            for (int k = 0; k < nd; ++k) x[k] = xdg[j + npe * k + npe * ncx * elem];
            for (int k = 0; k < nco; ++k) odg_local[k] = 0.0;

            PdeModel::cpuinitodg(odg_local, x, uinf, param);

            const int ncopy = (nco_runtime < nco) ? nco_runtime : nco;
            for (int k = 0; k < ncopy; ++k) {
                f[j + npe * k + npe * nco_runtime * elem] = odg_local[k];
            }
        }
    }
}
