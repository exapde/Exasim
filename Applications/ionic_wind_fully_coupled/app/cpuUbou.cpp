template <typename T> void cpuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param7 = param[6];
		T param9 = param[8];
		T param10 = param[9];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = -param7/(param9*param10);
	}
}

template <typename T> void cpuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
	}
}

template <typename T> void cpuUbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg8 = udg[7*ng+i];
		T udg12 = udg[11*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = t2+t3;
		T t5 = tanh(t4);
		T t6 = t5*1.0E+3;
		T t7 = t6-1.0;
		f[0*ng+i] = -t7*udg1;
		f[1*ng+i] = -t7*udg2;
		f[2*ng+i] = -t7*udg3;
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void cpuUbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg3 = udg[2*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = udg3;
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void cpuUbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg8 = udg[7*ng+i];
		T udg12 = udg[11*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = t2+t3;
		T t5 = tanh(t4);
		T t6 = t5*1.0E+3;
		T t7 = t6-1.0;
		f[0*ng+i] = -t7*udg1;
		f[1*ng+i] = -t7*udg2;
		f[2*ng+i] = -t7*udg3;
		f[3*ng+i] = udg4;
	}
}

template <typename T> void cpuUbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		cpuUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		cpuUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		cpuUbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		cpuUbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		cpuUbou5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void cpuUbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void cpuUbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
