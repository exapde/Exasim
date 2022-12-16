template <typename T> void opuUbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param4 = param[3];
		T param5 = param[4];
		T param7 = param[6];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T param14 = param[13];
		T uinf1 = uinf[0];
		T uinf3 = uinf[2];
		T xdg1 = xdg[0*ng+i];
		T t2 = 1.0/uinf1;
		T t3 = 1.0/uinf3;
		T t4 = xdg1*1.0E+2;
		T t5 = t4-5.0E+1;
		T t6 = (param1*t2)/2.0;
		T t7 = (param2*t2)/2.0;
		T t8 = (param3*t2)/2.0;
		T t9 = (param4*t2)/2.0;
		T t10 = (param5*t2)/2.0;
		T t11 = (param8*t2)/2.0;
		T t12 = (param7*t3)/2.0;
		T t13 = (param9*t2)/2.0;
		T t14 = (param10*t2)/2.0;
		T t15 = (param11*t2)/2.0;
		T t16 = (param12*t2)/2.0;
		T t17 = (param14*t3)/2.0;
		T t18 = tanh(t5);
		f[0*ng+i] = t6+t11-t18*(t6-t11);
		f[1*ng+i] = t7+t13-t18*(t7-t13);
		f[2*ng+i] = t8+t14-t18*(t8-t14);
		f[3*ng+i] = t9+t15-t18*(t9-t15);
		f[4*ng+i] = t10+t16-t18*(t10-t16);
		f[5*ng+i] = 0.0;
		f[6*ng+i] = t12+t17-t18*(t12-t17);
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
		T uinf1 = uinf[0];
		T uinf3 = uinf[2];
		T xdg1 = xdg[0*ng+i];
        T t1pi = 1.0/3.141592653589793;

		udg1 = udg1*(t1pi*atan(udg1*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg2 = udg2*(t1pi*atan(udg2*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg3 = udg3*(t1pi*atan(udg3*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg4 = udg4*(t1pi*atan(udg4*1.0E+3)+1.0/2.0)+3.183097800805168E-4;
		udg5 = udg5*(t1pi*atan(udg5*1.0E+3)+1.0/2.0)+3.183097800805168E-4;

		T Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7}; 

		f[0*ng+i] = Ucons[0];
		f[1*ng+i] = Ucons[1];
		f[2*ng+i] = Ucons[2];
		f[3*ng+i] = Ucons[3];
		f[4*ng+i] = Ucons[4];
		f[5*ng+i] = Ucons[5];
		uoutflow((double*) Ucons, (double) param[25], (double*) uinf, mix);

		f[6*ng+i] = Ucons[6];
	}
}

template <typename T> void opuUbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuUbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuUbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuUbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuUbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
