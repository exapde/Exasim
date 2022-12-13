template <typename T> void cpuTdfunc(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
}

template void cpuTdfunc(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuTdfunc(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
