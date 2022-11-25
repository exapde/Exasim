template <typename T> void opuSource2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T wdg1 = wdg[0*ng+i];
		f[0*ng+i] = wdg1;
	}
}

template void opuSource2(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuSource2(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
