template <typename T> void opuInitwdg(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne, Mutation::Mixture *mix)
{
}

template void opuInitwdg(double *, double *, double *, double *, int, int, int, int, int, int, Mutation::Mixture *);
template void opuInitwdg(float *, float *, float *, float *, int, int, int, int, int, int, Mutation::Mixture *);
