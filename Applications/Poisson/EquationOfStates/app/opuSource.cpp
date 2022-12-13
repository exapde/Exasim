template <typename T> void opuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		f[0*ng+i] = xdg1*sin(xdg2);
	}
}

template void opuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
