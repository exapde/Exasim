template <typename T> void cpuInitu(T *f, T *xdg, T *uinf, T *param, int modelnumber, int ng, int ncx, int nce, int npe, int ne)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param5 = param[4];
		T param7 = param[6];
		T param9 = param[8];
		T param10 = param[9];
		T param11 = param[10];
		T param12 = param[11];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		T t2 = xdg1*xdg1;
		T t3 = xdg2*xdg2;
		T t4 = t2+t3;
		T t5 = 1.0/param11;
		T t6 = param12*(1.0/4.0E1);
		T t7 = exp(t6);
		T t8 = param10-param11;
		T t9 = t7*t8;
		T t10 = sqrt(t4);
		T t11 = t10*(1.0/4.0E1);
		T t12 = 1.0/param1;
		T t17 = param7*t12;
		T t13 = param7-t17;
		T t14 = 1.0/t13;
		T t15 = param12*param12;
		T t16 = 1.0/t4;
		T t18 = exp(t11);
		T t19 = param11*t18;
		T t20 = t9+t19;
		T t21 = log(t20);
		T t22 = t5*t21*4.0E1;
		T t23 = param11*t7;
		T t24 = t9+t23;
		T t25 = log(t24);
		T t26 = t22-t5*t25*4.0E1;
		T t27 = exp(-param5*t14*t15*t16*t26);
		f[j+npe*0+npe*nce*k] = (param9*t14*t27)/(param11+t7*t8*exp(-t11));
		f[j+npe*1+npe*nce*k] = 0.0;
		f[j+npe*2+npe*nce*k] = 0.0;
		f[j+npe*3+npe*nce*k] = (param9*t27)/(param1-1.0);
	}
}

template void cpuInitu(double *, double *, double *, double *, int, int, int, int, int, int);
template void cpuInitu(float *, float *, float *, float *, int, int, int, int, int, int);
