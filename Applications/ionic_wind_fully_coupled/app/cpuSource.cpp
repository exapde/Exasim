template <typename T> void cpuSource(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param5 = param[4];
		T param6 = param[5];
		T param10 = param[9];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg8 = udg[7*ng+i];
		T udg12 = udg[11*ng+i];
		T t2 = udg8*udg8;
		T t3 = udg12*udg12;
		T t4 = 1.0/param5;
		T t5 = t2+t3;
		T t6 = sqrt(t5);
		f[0*ng+i] = param10*t6*udg1*(-2.2681E-19)+param6*t4*udg1*udg2*5.291005291005291E-12;
		f[1*ng+i] = param10*t6*udg1*1.19E-21+param6*t4*udg2*(udg1+udg3)*5.291005291005291E-12;
		f[2*ng+i] = param10*t6*udg1*2.28E-19+param6*t4*udg2*udg3*5.291005291005291E-12;
		f[3*ng+i] = -udg1+udg2-udg3;
	}
}

template void cpuSource(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void cpuSource(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
