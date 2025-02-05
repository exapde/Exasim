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
		dstype t16 = t2*2.0;
		dstype t17 = t3*2.0;
		dstype t4 = -t16-t17+2.0;
		dstype t5 = exp(t4);
		dstype t6 = param1*1.583143494411528E-1;
		dstype t7 = t6-1.583143494411528E-1;
		dstype t18 = t5*t7;
		dstype t8 = -t18+1.0;
		dstype t9 = param1-1.0;
		dstype t10 = 1.0/t9;
		dstype t11 = pow(t8,t10);
		dstype t12 = 1.0/param2;
		dstype t13 = 1.0/3.141592653589793;
		dstype t14 = -t2-t3+1.0;
		dstype t15 = exp(t14);
		dstype t19 = t12*t13*t15*xdg2*(5.0/2.0);
		dstype t20 = t19-1.0;
		dstype t21 = t10*2.0;
		dstype t22 = pow(t8,t21);
		dstype t23 = 1.0/(param2*param2);
		f[j+npe*0+npe*nce*k] = t11;
		f[j+npe*1+npe*nce*k] = -t11*t20;
		f[j+npe*2+npe*nce*k] = t11*t12*t13*t15*xdg1*(5.0/2.0);
		f[j+npe*3+npe*nce*k] = pow(t8,-t10)*((t20*t20)*t22*(1.0/2.0)+t2*t5*t22*t23*1.0/(3.141592653589793*3.141592653589793)*(2.5E1/8.0))+(t10*pow(t11,param1)*t23)/param1;
	}
}

