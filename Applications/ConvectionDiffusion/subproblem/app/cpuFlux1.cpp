template <typename T> void cpuFlux1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T wdg2 = wdg[1*ng+i];
		T wdg3 = wdg[2*ng+i];
		f[0*ng+i] = param1*udg2+udg1*wdg2;
		f[1*ng+i] = param1*udg3+udg1*wdg3;
	}
}

template void cpuFlux1(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuFlux1(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
