template <typename T> void gpuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int);
