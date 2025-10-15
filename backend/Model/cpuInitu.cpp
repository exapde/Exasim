void cpuInitu(dstype* f, const dstype* xdg, const dstype* uinf, const dstype* param, const int modelnumber, const int ng, const int ncx, const int nce, const int npe, const int ne)
{
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = i/npe;
		dstype param1 = param[0];
		dstype param2 = param[1];
		dstype xdg1 = xdg[j+npe*0+npe*ncx*k];
		dstype xdg2 = xdg[j+npe*1+npe*ncx*k];
		dstype t2 = xdg1*xdg1;
		dstype t3 = xdg2*xdg2;
		dstype t4 = 1.0/3.141592653589793;
		dstype t5 = param1-1.0;
		dstype t6 = 1.0/param2;
		dstype t22 = param1*1.583143494411528E-1;
		dstype t7 = t6*t6;
		dstype t8 = t2*2.0;
		dstype t9 = t3*2.0;
		dstype t10 = -t2;
		dstype t12 = -t3;
		dstype t14 = 1.0/t5;
		dstype t23 = t22-1.583143494411528E-1;
		dstype t11 = -t8;
		dstype t13 = -t9;
		dstype t15 = t14*2.0;
		dstype t16 = t10+t12+1.0;
		dstype t17 = t11+t13+2.0;
		dstype t18 = exp(t16);
		dstype t19 = exp(t17);
		dstype t20 = t4*t6*t18*xdg2*(5.0/2.0);
		dstype t21 = t20-1.0;
		dstype t24 = t19*t23;
		dstype t25 = -t24;
		dstype t26 = t25+1.0;
		dstype t27 = pow(t26,t14);
		dstype t28 = pow(t26,t15);
		f[j+npe*0+npe*nce*k] = t27;
		f[j+npe*1+npe*nce*k] = -t21*t27;
		f[j+npe*2+npe*nce*k] = t4*t6*t18*t27*xdg1*(5.0/2.0);
		f[j+npe*3+npe*nce*k] = (((t21*t21)*t28)/2.0+t2*(t4*t4)*t7*t19*t28*(2.5E+1/8.0))/t27+(t7*t14*pow(t27,param1))/param1;
	}
}

