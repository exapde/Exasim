template <typename T> void gpuInitq(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitq(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitq(float *, float *, float *, float *, int, int, int, int, int, int);
