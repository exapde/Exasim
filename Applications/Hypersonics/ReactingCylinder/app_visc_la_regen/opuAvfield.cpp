template <typename T> void opuAvfield(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
{
	double Ustate[6];
	int nspecies = 5;
	int ndim = 2;
	double rho_scale = uinf[0];
	double u_scale = uinf[1];
	double rhoe_scale = uinf[2];
	double omega_scale = rho_scale*u_scale;
    double mu_scale = uinf[4];
    double kappa_scale = uinf[5];
    double Uwork[5];
    double dTdU[6];
	double D_i[5];
    double h_i[5];
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param2 = param[1];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T param18 = param[17];
		T param22 = param[21];
		T param23 = param[22];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T xdg1 = xdg[j+npe*0+npe*ncx*k];
		T xdg2 = xdg[j+npe*1+npe*ncx*k];
		T udg1 = udg[j+npe*0+npe*nc*k];
		T udg2 = udg[j+npe*1+npe*nc*k];
		T udg3 = udg[j+npe*2+npe*nc*k];
		T udg4 = udg[j+npe*3+npe*nc*k];
		T udg5 = udg[j+npe*4+npe*nc*k];
		T udg6 = udg[j+npe*5+npe*nc*k];
		T udg7 = udg[j+npe*6+npe*nc*k];
		T udg8 = udg[j+npe*7+npe*nc*k];
		T udg9 = udg[j+npe*8+npe*nc*k];
		T udg10 = udg[j+npe*9+npe*nc*k];
		T udg11 = udg[j+npe*10+npe*nc*k];
		T udg12 = udg[j+npe*11+npe*nc*k];
		T udg13 = udg[j+npe*12+npe*nc*k];
		T udg14 = udg[j+npe*13+npe*nc*k];
		T udg15 = udg[j+npe*14+npe*nc*k];
		T udg16 = udg[j+npe*15+npe*nc*k];
		T udg17 = udg[j+npe*16+npe*nc*k];
		T udg18 = udg[j+npe*17+npe*nc*k];
		T udg19 = udg[j+npe*18+npe*nc*k];
		T udg20 = udg[j+npe*19+npe*nc*k];
		T udg21 = udg[j+npe*20+npe*nc*k];
		T udg22 = udg[j+npe*21+npe*nc*k];
		T udg23 = udg[j+npe*22+npe*nc*k];
		T udg24 = udg[j+npe*23+npe*nc*k];
		T odg1 = odg[j+npe*0+npe*nco*k];
		T odg2 = odg[j+npe*1+npe*nco*k];

		T t1pi = 1.0/3.141592653589793;
        
        udg1 = udg1*(t1pi*atan(udg1*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg2 = udg2*(t1pi*atan(udg2*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg3 = udg3*(t1pi*atan(udg3*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg4 = udg4*(t1pi*atan(udg4*1.0E+12)+1.0/2.0)+3.182454300088011E-13;
		udg5 = udg5*(t1pi*atan(udg5*1.0E+12)+1.0/2.0)+3.182454300088011E-13;

		double Ucons[8] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7, udg8};
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, ndim);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies, ndim);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

        mix->setState(rhovec, &rhoe, 0);

        uinf1 = mix->P() / rhoe_scale;
		uinf2 = mix->mixtureFrozenGamma();

		T t2 = atan(param22);
		T t3 = param22*udg1;
		T t4 = param22*udg2;
		T t5 = param22*udg3;
		T t6 = param22*udg4;
		T t7 = param22*udg5;
		T t8 = param23*2.0;
		T t9 = uinf2*2.0;
		T t10 = uinf2+1.0;
		T t16 = 1.0/3.141592653589793;
		T t17 = -param23;
		T t18 = 1.0/param15;
		T t21 = udg9+udg10+udg11+udg12+udg13;
		T t22 = udg17+udg18+udg19+udg20+udg21;
		T t11 = atan(t3);
		T t12 = atan(t4);
		T t13 = atan(t5);
		T t14 = atan(t6);
		T t15 = atan(t7);
		T t19 = t9-2.0;
		T t20 = 1.0/t10;
		T t23 = t2*t16;
		T t24 = t23*5.0;
		T t25 = t11*t16;
		T t26 = t12*t16;
		T t27 = t13*t16;
		T t28 = t14*t16;
		T t29 = t15*t16;
		T t30 = -t23;
		T t31 = -t24;
		T t32 = t25+1.0/2.0;
		T t33 = t26+1.0/2.0;
		T t34 = t27+1.0/2.0;
		T t35 = t28+1.0/2.0;
		T t36 = t29+1.0/2.0;
		T t37 = t32*udg1;
		T t38 = t33*udg2;
		T t39 = t34*udg3;
		T t40 = t35*udg4;
		T t41 = t36*udg5;
		T t42 = t31+t37+t38+t39+t40+t41+5.0/2.0;
		T t43 = 1.0/t42;
		T t44 = t43*t43;
		T t45 = t43*udg8;
		T t46 = t43*uinf1;
		T t47 = t21*t43*udg6;
		T t48 = t22*t43*udg7;
		T t49 = -t47;
		T t50 = -t48;
		T t53 = -t43*(t47-udg14);
		T t54 = -t43*(t48-udg23);
		T t55 = t45+t46-1.0E-4;
		T t61 = -param22*(-t17+t43*(t47-udg14)+t43*(t48-udg23));
		T t51 = t49+udg14;
		T t52 = t50+udg23;
		T t56 = param22*t55;
		T t60 = t17+t53+t54;
		T t62 = atan(t61);
		T t57 = atan(t56);
		T t63 = t16*t62;
		T t58 = t16*t57;
		T t65 = t63-1.0/2.0;
		T t59 = t58+1.0/2.0;
		T t69 = -t65*(-t17+t43*(t47-udg14)+t43*(t48-udg23));
		T t70 = t65*(-t17+t43*(t47-udg14)+t43*(t48-udg23));
		T t64 = t55*t59;
		T t71 = t8+t23+t70-1.0/2.0;
		T t66 = t30+t64+5.001E-1;
		T t72 = param22*t71;
		T t67 = t19*t20*t66;
		T t73 = atan(t72);
		T t68 = 1.0/sqrt(t67);
		T t74 = t16*t73;
		T t75 = t74-1.0/2.0;
		T t76 = t71*t75;
		T t77 = param23+t70+t76;
		T t78 = odg2*t18*t68*t77;
		T t79 = -t78;
		T t80 = param17+t79;
		T t81 = param22*t80;
		T t82 = atan(t81);
		T t83 = t16*t82;
		T t84 = t83-1.0/2.0;
		T t85 = t80*t84;
		T t86 = -t85;
		T t87 = param18+t23+t86-1.0/2.0;
		f[j+npe*0+npe*nce*k] = odg2*param16*t18*(sqrt(t44*(udg6*udg6)+t44*(udg7*udg7))+1.0/t68)*(t85-t87*(t16*atan(param22*t87)-1.0/2.0));
	}
}

template void opuAvfield(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuAvfield(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);