template <typename T> void gpuInitodg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
}

template void gpuInitodg(double *, double *, double *, double *, int, int, int, int, int, int);
template void gpuInitodg(float *, float *, float *, float *, int, int, int, int, int, int);
