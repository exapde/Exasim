template <typename T> void cpuFlux2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		f[0*ng+i] = udg2;
		f[1*ng+i] = udg3;
	}
}

template void cpuFlux2(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuFlux2(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
