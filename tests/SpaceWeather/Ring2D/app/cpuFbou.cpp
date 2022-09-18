template <typename T> void cpuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param3 = param[2];
		T param7 = param[6];
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
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = 1.0/udg1;
		T t3 = 1.0/(udg1*udg1);
		T t4 = udg2*udg2;
		T t5 = t3*t4*(1.0/2.0);
		T t6 = udg3*udg3;
		T t7 = t3*t6*(1.0/2.0);
		T t8 = t5+t7;
		T t11 = t8*udg1;
		T t9 = -t11+udg4;
		T t10 = param1-1.0;
		T t12 = 1.0/param1;
		T t23 = param7*t12;
		T t13 = param7-t23;
		T t14 = 1.0/(t13*t13);
		T t15 = t2*t9*t10*t14;
		T t16 = sqrt(t15);
		T t17 = t2*udg2*udg3;
		T t29 = t2*udg3*udg5;
		T t18 = -t29+udg7;
		T t19 = t2*t18;
		T t30 = t2*udg2*udg9;
		T t20 = -t30+udg10;
		T t21 = t2*t20;
		T t22 = t19+t21;
		T t24 = param2*t16*t22;
		T t25 = t17+t24;
		T t26 = t9*t10;
		T t31 = t2*udg2*udg5;
		T t27 = -t31+udg6;
		T t33 = t2*udg3*udg9;
		T t28 = -t33+udg11;
		T t32 = t2*t27*2.0;
		T t34 = t32-t2*t28;
		T t35 = 1.0/t13;
		T t36 = t2*udg4;
		T t37 = t2*t9*t10;
		T t38 = t36+t37;
		T t39 = t2*t27;
		T t40 = t39-t2*t28*2.0;
		T t41 = t2*t9*t10*t35;
		T t42 = pow(t41,3.0/4.0);
		f[0*ng+i] = 0.0;
		f[1*ng+i] = nlg2*t25+tau1*(udg2-uhg2)+nlg1*(t26+t2*t4+param2*t16*t34*(2.0/3.0));
		f[2*ng+i] = nlg1*t25+tau1*(udg3-uhg3)+nlg2*(t26+t2*t6-param2*t16*t40*(2.0/3.0));
		f[3*ng+i] = nlg1*(t38*udg2+param2*t2*t16*t22*udg3+param2*t2*t16*t34*udg2*(2.0/3.0)-param3*t3*t35*t42*(t10*udg1*(-udg8+t8*udg5+udg1*(t3*t18*udg3+t3*t27*udg2))+t9*t10*udg5))+nlg2*(t38*udg3+param2*t2*t16*t22*udg2-param2*t2*t16*t40*udg3*(2.0/3.0)-param3*t3*t35*t42*(t10*udg1*(-udg12+t8*udg9+udg1*(t3*t20*udg2+t3*t28*udg3))+t9*t10*udg9))+tau1*(udg4-uhg4);
	}
}

template <typename T> void cpuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param2 = param[1];
		T param7 = param[6];
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
		T t3 = 1.0/(udg1*udg1);
		T t4 = udg2*udg2;
		T t5 = t3*t4*(1.0/2.0);
		T t6 = udg3*udg3;
		T t7 = t3*t6*(1.0/2.0);
		T t8 = t5+t7;
		T t11 = t8*udg1;
		T t9 = -t11+udg4;
		T t10 = param1-1.0;
		T t12 = 1.0/param1;
		T t23 = param7*t12;
		T t13 = param7-t23;
		T t14 = 1.0/(t13*t13);
		T t15 = t2*t9*t10*t14;
		T t16 = sqrt(t15);
		T t17 = t2*udg2*udg3;
		T t18 = udg7-t2*udg3*udg5;
		T t19 = t2*t18;
		T t20 = udg10-t2*udg2*udg9;
		T t21 = t2*t20;
		T t22 = t19+t21;
		T t24 = param2*t16*t22;
		T t25 = t17+t24;
		T t26 = t9*t10;
		T t27 = udg6-t2*udg2*udg5;
		T t28 = udg11-t2*udg3*udg9;
		f[0*ng+i] = 0.0;
		f[1*ng+i] = nlg1*(t26+t2*t4+param2*t16*(t2*t27*2.0-t2*t28)*(2.0/3.0))+nlg2*t25+tau1*(udg2-uhg2);
		f[2*ng+i] = nlg2*(t26+t2*t6-param2*t16*(t2*t27-t2*t28*2.0)*(2.0/3.0))+nlg1*t25+tau1*(udg3-uhg3);
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void cpuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		cpuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		cpuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void cpuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void cpuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
