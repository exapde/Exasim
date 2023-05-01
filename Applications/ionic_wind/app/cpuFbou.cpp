template <typename T> void cpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		f[0*ng+i] = 0.0;
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*xdg1*(odg2+udg4)+nlg2*xdg1*(odg3+udg6);
	}
}

template <typename T> void cpuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
	}
}

template <typename T> void cpuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = odg2+udg4;
		T t3 = odg3+udg6;
		T t4 = -uhg1;
		T t5 = t4+udg1;
		T t7 = nlg1*t2*1.0E+3;
		T t8 = nlg2*t3*1.0E+3;
		T t6 = t5*tau1;
		T t9 = t7+t8;
		T t10 = tanh(t9);
		T t11 = t10/2.0;
		f[0*ng+i] = -(t6-udg1*xdg1*(nlg1*t2+nlg2*t3))*(t11-1.0/2.0)+(t11+1.0/2.0)*(t6+nlg1*xdg1*(udg3/1.0E+1-t2*udg1)+nlg2*xdg1*(udg5/1.0E+1-t3*udg1));
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t2*xdg1+nlg2*t3*xdg1;
	}
}

template <typename T> void cpuFbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T odg2 = odg[1*ng+i];
		T odg3 = odg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = odg2+udg4;
		T t3 = odg3+udg6;
		f[0*ng+i] = tau1*(udg1-uhg1)-udg1*xdg1*(nlg1*t2+nlg2*t3);
		f[1*ng+i] = tau1*(udg2-uhg2)+nlg1*t2*xdg1+nlg2*t3*xdg1;
	}
}

template <typename T> void cpuFbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
	}
}

template <typename T> void cpuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		cpuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		cpuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		cpuFbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		cpuFbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		cpuFbou5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void cpuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void cpuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
