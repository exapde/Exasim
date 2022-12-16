template <typename T> void opuInitodg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne, Mutation::Mixture *mix)
{
}

template void opuInitodg(double *, double *, double *, double *, int, int, int, int, int, int, Mutation::Mixture *);
template void opuInitodg(float *, float *, float *, float *, int, int, int, int, int, int, Mutation::Mixture *);
