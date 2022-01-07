template <typename T> void cpuTdfunc(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param13 = param[12];
		T param14 = param[13];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T xdg1 = xdg[0*ng+i];
		T xdg2 = xdg[1*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T odg1 = odg[0*ng+i];
		T odg2 = odg[1*ng+i];
		f[0*ng+i] = 1.0;
		f[1*ng+i] = 1.0;
		f[2*ng+i] = 1.0;
		f[3*ng+i] = 1.0;
	}
}

template void cpuTdfunc(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int);
template void cpuTdfunc(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int);
