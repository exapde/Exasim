template <typename T> void gpuInitq1(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitq1(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitq1(float *, float *, float *, float *, int, int, int, int, int, int);
