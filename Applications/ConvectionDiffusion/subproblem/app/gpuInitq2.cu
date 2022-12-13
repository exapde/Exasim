template <typename T> void gpuInitq2(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitq2(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitq2(float *, float *, float *, float *, int, int, int, int, int, int);
