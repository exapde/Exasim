template <typename T> void gpuInitodg1(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitodg1(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitodg1(float *, float *, float *, float *, int, int, int, int, int, int);
