template <typename T> void opuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T t2 = odg2+udg4;
		T t3 = odg3+udg6;
		f[0*ng+i] = xdg1*(udg3/1.0E+1-t2*udg1);
		f[1*ng+i] = t2*xdg1;
		f[2*ng+i] = xdg1*(udg5/1.0E+1-t3*udg1);
		f[3*ng+i] = t3*xdg1;
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int);
