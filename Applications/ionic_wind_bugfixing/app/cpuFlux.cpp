template <typename T> void cpuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param5 = param[4];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T t2 = odg2+udg8;
		T t3 = odg3+udg12;
		f[0*ng+i] = xdg1*(param5*udg5+t2*udg1);
		f[1*ng+i] = xdg1*(param5*udg6+t2*udg2);
		f[2*ng+i] = xdg1*(param5*udg7+t2*udg3);
		f[3*ng+i] = t2*xdg1;
		f[4*ng+i] = xdg1*(param5*udg9+t3*udg1);
		f[5*ng+i] = xdg1*(param5*udg10+t3*udg2);
		f[6*ng+i] = xdg1*(param5*udg11+t3*udg3);
		f[7*ng+i] = t3*xdg1;
	}
}

template void cpuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);