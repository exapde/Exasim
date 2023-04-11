template <typename T> void opuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
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
		T t2 = 1.0/param9;
		T t3 = 1.0/param10;
		f[0*ng+i] = -xdg1*(udg1*udg8+t2*t3*udg5*(8.0E+2/1.89E+2));
		f[1*ng+i] = -xdg1*(udg2*udg8*6.19047619047619E-3+t2*t3*udg6*(8.0E+2/1.89E+2));
		f[2*ng+i] = -xdg1*((udg3*udg8)/1.4E+2+t2*t3*udg7*(8.0E+2/1.89E+2));
		f[3*ng+i] = udg8*xdg1;
		f[4*ng+i] = -xdg1*(udg1*udg12+t2*t3*udg9*(8.0E+2/1.89E+2));
		f[5*ng+i] = -xdg1*(udg2*udg12*6.19047619047619E-3+t2*t3*udg10*(8.0E+2/1.89E+2));
		f[6*ng+i] = -xdg1*((udg3*udg12)/1.4E+2+t2*t3*udg11*(8.0E+2/1.89E+2));
		f[7*ng+i] = udg12*xdg1;
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
