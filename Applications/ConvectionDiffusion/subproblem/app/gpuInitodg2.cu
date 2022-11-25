template <typename T> void gpuInitodg2(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitodg2(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitodg2(float *, float *, float *, float *, int, int, int, int, int, int);
