template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param8 = param[7];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg7 = udg[6*ng+i];
		T udg8 = udg[7*ng+i];
		T udg11 = udg[10*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = t2+t3;
		f[0*ng+i] = param9*param10*(t6*udg1-param8*udg2*sqrt(udg8*udg8+udg12*udg12))*(1.89E+2/8.0E+2);
		f[1*ng+i] = t6*udg2*(-6.19047619047619E-3)+tau1*(udg2-uhg2);
		f[2*ng+i] = tau1*(udg3-uhg3)-nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))-nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2));
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
	}
}

template <typename T> void opuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void opuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
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
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = -uhg1;
		T t7 = -uhg2;
		T t8 = -uhg3;
		T t9 = t6+udg1;
		T t10 = t7+udg2;
		T t11 = t8+udg3;
		T t12 = t2+t3;
		T t13 = tanh(t12);
		T t14 = t9*tau1;
		T t15 = t10*tau1;
		T t16 = t11*tau1;
		T t17 = t13*1.0E+3;
		T t18 = t17-1.0;
		f[0*ng+i] = t13*(-t14+nlg1*xdg1*(udg1*udg8+t4*t5*udg5*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg1*udg12+t4*t5*udg9*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t14-t12*udg1);
		f[1*ng+i] = t13*(-t15+nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t15-t12*udg2*6.19047619047619E-3);
		f[2*ng+i] = t13*(-t16+nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))+nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t16-(t12*udg3)/1.4E+2);
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
	}
}

template <typename T> void opuFbou4(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
		T xdg1 = xdg[0*ng+i];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg6 = udg[5*ng+i];
		T udg8 = udg[7*ng+i];
		T udg10 = udg[9*ng+i];
		T udg12 = udg[11*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = t2+t3;
		f[0*ng+i] = -t6*udg1+tau1*(udg1-uhg1);
		f[1*ng+i] = tau1*(udg2-uhg2)-nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))-nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2));
		f[2*ng+i] = t6*udg3*(-1.0/1.4E+2)+tau1*(udg3-uhg3);
		f[3*ng+i] = t2*xdg1+t3*xdg1+tau1*(udg4-uhg4);
	}
}

template <typename T> void opuFbou5(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	for (int i = 0; i <ng; i++) {
		T param1 = param[0];
		T param9 = param[8];
		T param10 = param[9];
		T tau1 = tau[0];
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
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
		T t2 = nlg1*udg8;
		T t3 = nlg2*udg12;
		T t4 = 1.0/param9;
		T t5 = 1.0/param10;
		T t6 = -uhg1;
		T t7 = -uhg2;
		T t8 = -uhg3;
		T t9 = t6+udg1;
		T t10 = t7+udg2;
		T t11 = t8+udg3;
		T t12 = t2+t3;
		T t13 = tanh(t12);
		T t14 = t9*tau1;
		T t15 = t10*tau1;
		T t16 = t11*tau1;
		T t17 = t13*1.0E+3;
		T t18 = t17-1.0;
		f[0*ng+i] = t13*(-t14+nlg1*xdg1*(udg1*udg8+t4*t5*udg5*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg1*udg12+t4*t5*udg9*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t14-t12*udg1);
		f[1*ng+i] = t13*(-t15+nlg1*xdg1*(udg2*udg8*6.19047619047619E-3+t4*t5*udg6*(8.0E+2/1.89E+2))+nlg2*xdg1*(udg2*udg12*6.19047619047619E-3+t4*t5*udg10*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t15-t12*udg2*6.19047619047619E-3);
		f[2*ng+i] = t13*(-t16+nlg1*xdg1*((udg3*udg8)/1.4E+2+t4*t5*udg7*(8.0E+2/1.89E+2))+nlg2*xdg1*((udg3*udg12)/1.4E+2+t4*t5*udg11*(8.0E+2/1.89E+2)))*-1.0E+3-t18*(t16-(t12*udg3)/1.4E+2);
		f[3*ng+i] = 0.0;
	}
}

template <typename T> void opuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw)
{
	if (ib == 1)
		opuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 2)
		opuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 3)
		opuFbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 4)
		opuFbou4(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
	else if (ib == 5)
		opuFbou5(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw);
}

template void opuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int);
template void opuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int);
