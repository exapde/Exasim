template <typename T> void opuAvfield(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, int nce, int npe, int ne, Mutation::Mixture *mix)
{
    double Ustate[6];
	int nspecies = 5;
	int ndim = 2;
	double rho_scale = uinf[0];
	double u_scale = uinf[1];
	double rhoe_scale = uinf[2];
	double omega_scale = rho_scale*u_scale;
	for (int i = 0; i <ng; i++) {
		int j = i%npe;
		int k = (i-j)/npe;
		T param1 = param[0];
		T param15 = param[14];
		T param16 = param[15];
		T param17 = param[16];
		T param18 = param[17];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf3 = uinf[2];
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
		T odg3 = odg[j+npe*2+npe*nco*k];
        
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
        // uinf3 = (mix->frozenThermalConductivity() / uinf[5]) / (mix->viscosity() / uinf[4]); // cp / Pr^*
		uinf3 = (mix->mixtureFrozenCpMass() / 1182.1920097928833) / 0.9;

		T t2 = udg6*udg6;
		T t3 = udg7*udg7;
		T t4 = uinf2*2.0;
		T t5 = uinf2+1.0;
		T t6 = 1.0/3.141592653589793;
		T t7 = -param18;
		T t8 = 1.0/param15;
		T t11 = udg9+udg10+udg11+udg12+udg13;
		T t12 = udg17+udg18+udg19+udg20+udg21;
		T t13 = param17*1.0E+12;
		T t14 = param18*1.0E+12;
		T t15 = udg1*1.0E+12;
		T t16 = udg2*1.0E+12;
		T t17 = udg3*1.0E+12;
		T t18 = udg4*1.0E+12;
		T t19 = udg5*1.0E+12;
		T t9 = t4-2.0;
		T t10 = 1.0/t5;
		T t20 = atan(t15);
		T t21 = atan(t16);
		T t22 = atan(t17);
		T t23 = atan(t18);
		T t24 = atan(t19);
		T t25 = -t14;
		T t26 = t6*t20;
		T t27 = t6*t21;
		T t28 = t6*t22;
		T t29 = t6*t23;
		T t30 = t6*t24;
		T t31 = t26+1.0/2.0;
		T t32 = t27+1.0/2.0;
		T t33 = t28+1.0/2.0;
		T t34 = t29+1.0/2.0;
		T t35 = t30+1.0/2.0;
		T t36 = t31*udg1;
		T t37 = t32*udg2;
		T t38 = t33*udg3;
		T t39 = t34*udg4;
		T t40 = t35*udg5;
		T t41 = t36+t37+t38+t39+t40+1.591227150044006E-12;
		T t42 = 1.0/t41;
		T t43 = t42*t42;
		T t44 = t42*udg8;
		T t45 = t42*uinf1;
		T t48 = t11*t42*udg6;
		T t49 = t12*t42*udg7;
		T t46 = t2*t43;
		T t47 = t3*t43;
		T t50 = -t48;
		T t51 = -t49;
		T t54 = t44*1.0E+12;
		T t55 = t45*1.0E+12;
		T t56 = -t42*(t48-udg14);
		T t57 = -t42*(t49-udg23);
		T t58 = t44+t45-1.0E-4;
		T t59 = t42*(t48-udg14)*-1.0E+12;
		T t60 = t42*(t49-udg23)*-1.0E+12;
		T t52 = t50+udg14;
		T t53 = t51+udg23;
		T t61 = t54+t55-1.0E+8;
		T t65 = t56+t57-1.0E+4;
		T t67 = t59+t60-1.0E+16;
		T t62 = atan(t61);
		T t68 = atan(t67);
		T t63 = t6*t62;
		T t69 = t6*t68;
		T t64 = t63+1.0/2.0;
		T t70 = t69-1.0/2.0;
		T t66 = t58*t64;
		T t76 = -t70*(t42*(t48-udg14)+t42*(t49-udg23)+1.0E+4);
		T t77 = t70*(t42*(t48-udg14)+t42*(t49-udg23)+1.0E+4);
		T t71 = t66+1.000000003182454E-4;
		T t78 = t77*-1.0E+12;
		T t79 = t76-2.0E+4;
		T t72 = t9*t10*t71;
		T t80 = t78-2.0E+16;
		T t73 = 1.0/sqrt(t72);
		T t74 = t46+t47+t72;
		T t81 = atan(t80);
		T t75 = sqrt(t74);
		T t82 = t6*t81;
		T t83 = t82+1.0/2.0;
		T t84 = -t83*(t77+2.0E+4);
		T t85 = t77+t84+1.0E+4;
		T t86 = odg3*t8*t73*t85;
		T t87 = -t86;
		T t89 = t86*1.0E+12;
		T t88 = param17+t87;
		T t90 = -t89;
		T t91 = t13+t90;
		T t92 = atan(t91);
		T t93 = t6*t92;
		T t94 = t93-1.0/2.0;
		T t95 = t88*t94;
		T t96 = t95*1.0E+12;
		T t97 = t7+t95+3.182454300088011E-13;
		T t98 = t25+t96+3.182454300088011E-1;
		T t99 = atan(t98);
		T t100 = t6*t99;
		T t101 = t100+1.0/2.0;
		T t102 = t97*t101;
		T t103 = -t102;
		T t104 = t95+t103;
		f[j+npe*0+npe*nce*k] = odg3*param16*t8*t41*t75*t104;
		f[j+npe*1+npe*nce*k] = 0.5 * odg3*param16*t8*t41*t75*t104*uinf3;
	}
}

template void opuAvfield(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuAvfield(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);