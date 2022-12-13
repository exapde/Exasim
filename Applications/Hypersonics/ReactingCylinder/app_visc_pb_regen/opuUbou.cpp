template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param6 = param[5];
		T param7 = param[6];
		T param8 = param[7];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf3 = uinf[2];
		T t2 = 1.0/uinf1;
		T t3 = 1.0/uinf2;
		f[0*ng+i] = param1*t2;
		f[1*ng+i] = param2*t2;
		f[2*ng+i] = param3*t2;
		f[3*ng+i] = param4*t2;
		f[4*ng+i] = param5*t2;
		f[5*ng+i] = param6*t2*t3;
		f[6*ng+i] = param7*t2*t3;
		f[7*ng+i] = param8/uinf3;
	}
}

template <typename T> void opuUbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
		f[4*ng+i] = udg5;
		f[5*ng+i] = udg6;
		f[6*ng+i] = udg7;
		f[7*ng+i] = udg8;
	}
}

template <typename T> void opuUbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg8 = udg[7*ng+i];
		f[0*ng+i] = udg1;
		f[1*ng+i] = udg2;
		f[2*ng+i] = udg3;
		f[3*ng+i] = udg4;
		f[4*ng+i] = udg5;
		f[5*ng+i] = 0.0;
		f[6*ng+i] = 0.0;
		f[7*ng+i] = udg8;
	}
}

template <typename T> void opuUbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 3)
		opuUbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuUbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuUbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
