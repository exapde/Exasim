template <typename T> void gpuTdfunc1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
}

template void gpuTdfunc1(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void gpuTdfunc1(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);

template <typename T> void __device__ devicegpuTdfunc1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
}

