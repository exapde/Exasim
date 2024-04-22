template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T tau1 = tau[0];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		T udg9 = udg[8*ng+i];
		T udg10 = udg[9*ng+i];
		T udg11 = udg[10*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = 1.0/udg1;
		T t3 = udg2*udg2;
		T t4 = 1.0/(udg1*udg1);
		T t5 = 1.0/param2;
		T t6 = udg7-t2*udg3*udg5;
		T t7 = t2*t6;
		T t8 = udg10-t2*udg2*udg9;
		T t9 = t2*t8;
		T t10 = t7+t9;
		T t11 = t5*t10;
		T t12 = t2*udg2*udg3;
		T t13 = t11+t12;
		T t14 = udg3*udg3;
		T t15 = t3*t4*(1.0/2.0);
		T t16 = t4*t14*(1.0/2.0);
		T t17 = t15+t16;
		T t18 = udg4-t17*udg1;
		T t19 = param1-1.0;
		T t20 = t18*t19;
		T t21 = udg6-t2*udg2*udg5;
		T t22 = udg11-t2*udg3*udg9;
		f[0*ng+i] = 0.0;
		f[1*ng+i] = nlg2*t13+tau1*(udg2-uhg2)+nlg1*(t20+t2*t3+t5*(t2*t21*2.0-t2*t22)*(2.0/3.0));
		f[2*ng+i] = nlg1*t13+tau1*(udg3-uhg3)+nlg2*(t20+t2*t14-t5*(t2*t21-t2*t22*2.0)*(2.0/3.0));
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void opuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T tau1 = tau[0];
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
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = 1.0/udg1;
		T t3 = udg2*udg2;
		T t4 = 1.0/(udg1*udg1);
		T t5 = 1.0/param2;
		T t24 = t2*udg3*udg5;
		T t6 = -t24+udg7;
		T t7 = t2*t6;
		T t25 = t2*udg2*udg9;
		T t8 = -t25+udg10;
		T t9 = t2*t8;
		T t10 = t7+t9;
		T t11 = t5*t10;
		T t12 = t2*udg2*udg3;
		T t13 = t11+t12;
		T t14 = udg3*udg3;
		T t15 = t3*t4*(1.0/2.0);
		T t16 = t4*t14*(1.0/2.0);
		T t17 = t15+t16;
		T t23 = t17*udg1;
		T t18 = -t23+udg4;
		T t19 = param1-1.0;
		T t20 = t18*t19;
		T t26 = t2*udg2*udg5;
		T t21 = -t26+udg6;
		T t28 = t2*udg3*udg9;
		T t22 = -t28+udg11;
		T t27 = t2*t21*2.0;
		T t29 = t27-t2*t22;
		T t30 = t2*udg4;
		T t31 = t2*t18*t19;
		T t32 = t30+t31;
		T t33 = t2*t21;
		T t34 = t33-t2*t22*2.0;
		T t35 = 1.0/param3;
		T t36 = 1.0/t19;
		f[0*ng+i] = nlg1*udg2+nlg2*udg3+tau1*(udg1-uhg1);
		f[1*ng+i] = nlg1*(t20+t2*t3+t5*t29*(2.0/3.0))+nlg2*t13+tau1*(udg2-uhg2);
		f[2*ng+i] = nlg2*(t20+t2*t14-t5*t34*(2.0/3.0))+nlg1*t13+tau1*(udg3-uhg3);
		f[3*ng+i] = nlg1*(t32*udg2+t2*t5*t10*udg3+t2*t5*t29*udg2*(2.0/3.0)-param1*t4*t5*t35*t36*(t19*udg1*(-udg8+t17*udg5+udg1*(t4*t6*udg3+t4*t21*udg2))+t18*t19*udg5))+nlg2*(t32*udg3+t2*t5*t10*udg2-t2*t5*t34*udg3*(2.0/3.0)-param1*t4*t5*t35*t36*(t19*udg1*(-udg12+t17*udg9+udg1*(t4*t8*udg2+t4*t22*udg3))+t18*t19*udg9))+tau1*(udg4-uhg4);
	}
}

template <typename T> void opuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		opuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		opuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void opuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void opuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
