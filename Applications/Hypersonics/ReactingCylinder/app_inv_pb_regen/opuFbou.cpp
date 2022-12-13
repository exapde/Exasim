
template <typename T> void opuFbou1(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
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
		T param2 = param[1];
		T param21 = param[20];
		T param22 = param[21];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf14 = uinf[13];
		T uinf15 = uinf[14];
		T uinf16 = uinf[15];
		T uinf17 = uinf[16];
		T uinf18 = uinf[17];
		T uinf19 = uinf[18];
		T uinf20 = uinf[19];
		T tau1 = tau[0];
		T tau2 = tau[1];
		T tau3 = tau[2];
		T tau4 = tau[3];
		T tau5 = tau[4];
		T tau6 = tau[5];
		T tau7 = tau[6];
		T tau8 = tau[7];
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
		T udg13 = udg[12*ng+i];
		T udg14 = udg[13*ng+i];
		T udg15 = udg[14*ng+i];
		T udg16 = udg[15*ng+i];
		T udg17 = udg[16*ng+i];
		T udg18 = udg[17*ng+i];
		T udg19 = udg[18*ng+i];
		T udg20 = udg[19*ng+i];
		T udg21 = udg[20*ng+i];
		T udg22 = udg[21*ng+i];
		T udg23 = udg[22*ng+i];
		T udg24 = udg[23*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T uhg8 = uhg[7*ng+i];
		T odg1 = odg[0*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];


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

		uinf20 = 1.0/0.9;

        dT_dUstate(dTdU, Ustate, Uwork, nspecies, ndim, mix);
        nondimensionalize_dT_dUstate(dTdU, (double*)uinf, nspecies, nd);
        //uinf14 = 0.0; uinf15 = 0.0; uinf16 = 0.0; uinf17 = 0.0; uinf18 = 0.0; uinf19 = 0.0; // LAPLACIAN
		uinf14 = dTdU[0]; uinf15 = dTdU[1]; uinf16 = dTdU[2]; uinf17 = dTdU[3]; uinf18 = dTdU[4]; uinf19 = dTdU[5];

		T t2 = atan(param22);
		T t3 = param22*udg1;
		T t4 = param22*udg2;
		T t5 = param22*udg3;
		T t6 = param22*udg4;
		T t7 = param22*udg5;
		T t8 = param22*param22;
		T t9 = udg1*udg1;
		T t10 = udg2*udg2;
		T t11 = udg3*udg3;
		T t12 = udg4*udg4;
		T t13 = udg5*udg5;
		T t14 = udg6*udg6;
		T t15 = udg7*udg7;
		T t21 = 1.0/3.141592653589793;
		T t16 = atan(t3);
		T t17 = atan(t4);
		T t18 = atan(t5);
		T t19 = atan(t6);
		T t20 = atan(t7);
		T t22 = t3*t3;
		T t23 = t4*t4;
		T t24 = t5*t5;
		T t25 = t6*t6;
		T t26 = t7*t7;
		T t27 = t2*t21;
		T t28 = t27*5.0;
		T t29 = t22+1.0;
		T t30 = t23+1.0;
		T t31 = t24+1.0;
		T t32 = t25+1.0;
		T t33 = t26+1.0;
		T t34 = t16*t21;
		T t35 = t17*t21;
		T t36 = t18*t21;
		T t37 = t19*t21;
		T t38 = t20*t21;
		T t39 = -t27;
		T t40 = -t28;
		T t41 = 1.0/t29;
		T t42 = 1.0/t30;
		T t43 = 1.0/t31;
		T t44 = 1.0/t32;
		T t45 = 1.0/t33;
		T t46 = t34+1.0/2.0;
		T t47 = t35+1.0/2.0;
		T t48 = t36+1.0/2.0;
		T t49 = t37+1.0/2.0;
		T t50 = t38+1.0/2.0;
		T t51 = t46*udg1;
		T t52 = t47*udg2;
		T t53 = t48*udg3;
		T t54 = t49*udg4;
		T t55 = t50*udg5;
		T t56 = t3*t21*t41;
		T t57 = t4*t21*t42;
		T t58 = t5*t21*t43;
		T t59 = t6*t21*t44;
		T t60 = t7*t21*t45;
		T t61 = t39+t51+1.0/2.0;
		T t62 = t39+t52+1.0/2.0;
		T t63 = t39+t53+1.0/2.0;
		T t64 = t39+t54+1.0/2.0;
		T t65 = t39+t55+1.0/2.0;
		T t66 = t46+t56;
		T t67 = t47+t57;
		T t68 = t48+t58;
		T t69 = t49+t59;
		T t70 = t50+t60;
		T t81 = t40+t51+t52+t53+t54+t55+5.0/2.0;
		T t71 = t66*udg9;
		T t72 = t67*udg10;
		T t73 = t68*udg11;
		T t74 = t66*udg17;
		T t75 = t69*udg12;
		T t76 = t67*udg18;
		T t77 = t70*udg13;
		T t78 = t68*udg19;
		T t79 = t69*udg20;
		T t80 = t70*udg21;
		T t82 = 1.0/t81;
		T t83 = t82*t82;
		T t84 = t82*udg8;
		T t85 = t82*uinf1;
		T t88 = t71+t72+t73+t75+t77;
		T t89 = t74+t76+t78+t79+t80;
		T t86 = (t14*t83)/2.0;
		T t87 = (t15*t83)/2.0;
		T t90 = t84+t85;
		T t92 = t82*t88*udg6;
		T t93 = t82*t89*udg7;
		T t91 = t86+t87;
		T t94 = -t92;
		T t95 = -t93;
		T t98 = -t82*(t92-udg14);
		T t99 = -t82*(t93-udg23);
		T t101 = -odg1*(t82*(t92-udg14)+t82*(t93-udg23));
		T t96 = t94+udg14;
		T t97 = t95+udg23;
		T t100 = t98+t99;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t61*t82*udg6+nlg2*t61*t82*udg7;
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t62*t82*udg6+nlg2*t62*t82*udg7;
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t63*t82*udg6+nlg2*t63*t82*udg7;
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t64*t82*udg6+nlg2*t64*t82*udg7;
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t65*t82*udg6+nlg2*t65*t82*udg7;
		f[5*ng+i] = nlg1*(t101+uinf1+t14*t82)+tau6*(udg6-uhg6)+nlg2*t82*udg6*udg7;
		f[6*ng+i] = nlg2*(t101+uinf1+t15*t82)+tau7*(udg7-uhg7)+nlg1*t82*udg6*udg7;
		f[7*ng+i] = nlg1*(t90*udg6+t82*t101*udg6+odg1*param21*uinf20*(t71*uinf14+t72*uinf15+t73*uinf16+t75*uinf17+t77*uinf18-uinf19*(-udg16+t88*t91+t81*(t83*udg7*(udg15-t82*t88*udg7)-t83*udg6*(t92-udg14)))))+nlg2*(t90*udg7+t82*t101*udg7+odg1*param21*uinf20*(t74*uinf14+t76*uinf15+t78*uinf16+t79*uinf17+t80*uinf18-uinf19*(-udg24+t89*t91+t81*(t83*udg6*(udg22-t82*t89*udg6)-t83*udg7*(t93-udg23)))))+tau8*(udg8-uhg8);
	}
}

template <typename T> void opuFbou2(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
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
		T param2 = param[1];
		T param21 = param[20];
		T param22 = param[21];
		T uinf1 = uinf[0];
		T uinf2 = uinf[1];
		T uinf14 = uinf[13];
		T uinf15 = uinf[14];
		T uinf16 = uinf[15];
		T uinf17 = uinf[16];
		T uinf18 = uinf[17];
		T uinf19 = uinf[18];
		T uinf20 = uinf[19];
		T tau1 = tau[0];
		T tau2 = tau[1];
		T tau3 = tau[2];
		T tau4 = tau[3];
		T tau5 = tau[4];
		T tau6 = tau[5];
		T tau7 = tau[6];
		T tau8 = tau[7];
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
		T udg13 = udg[12*ng+i];
		T udg14 = udg[13*ng+i];
		T udg15 = udg[14*ng+i];
		T udg16 = udg[15*ng+i];
		T udg17 = udg[16*ng+i];
		T udg18 = udg[17*ng+i];
		T udg19 = udg[18*ng+i];
		T udg20 = udg[19*ng+i];
		T udg21 = udg[20*ng+i];
		T udg22 = udg[21*ng+i];
		T udg23 = udg[22*ng+i];
		T udg24 = udg[23*ng+i];
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T uhg8 = uhg[7*ng+i];
		T odg1 = odg[0*ng+i];
		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];

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

		uinf20 = 1.0/0.9;

        dT_dUstate(dTdU, Ustate, Uwork, nspecies, ndim, mix);
        nondimensionalize_dT_dUstate(dTdU, (double*)uinf, nspecies, nd);
        //uinf14 = 0.0; uinf15 = 0.0; uinf16 = 0.0; uinf17 = 0.0; uinf18 = 0.0; uinf19 = 0.0; // LAPLACIAN
		uinf14 = dTdU[0]; uinf15 = dTdU[1]; uinf16 = dTdU[2]; uinf17 = dTdU[3]; uinf18 = dTdU[4]; uinf19 = dTdU[5];

		T t2 = atan(param22);
		T t3 = param22*udg1;
		T t4 = param22*udg2;
		T t5 = param22*udg3;
		T t6 = param22*udg4;
		T t7 = param22*udg5;
		T t8 = param22*param22;
		T t9 = udg1*udg1;
		T t10 = udg2*udg2;
		T t11 = udg3*udg3;
		T t12 = udg4*udg4;
		T t13 = udg5*udg5;
		T t14 = udg6*udg6;
		T t15 = udg7*udg7;
		T t21 = 1.0/3.141592653589793;
		T t16 = atan(t3);
		T t17 = atan(t4);
		T t18 = atan(t5);
		T t19 = atan(t6);
		T t20 = atan(t7);
		T t22 = t3*t3;
		T t23 = t4*t4;
		T t24 = t5*t5;
		T t25 = t6*t6;
		T t26 = t7*t7;
		T t27 = t2*t21;
		T t28 = t27*5.0;
		T t29 = t22+1.0;
		T t30 = t23+1.0;
		T t31 = t24+1.0;
		T t32 = t25+1.0;
		T t33 = t26+1.0;
		T t34 = t16*t21;
		T t35 = t17*t21;
		T t36 = t18*t21;
		T t37 = t19*t21;
		T t38 = t20*t21;
		T t39 = -t27;
		T t40 = -t28;
		T t41 = 1.0/t29;
		T t42 = 1.0/t30;
		T t43 = 1.0/t31;
		T t44 = 1.0/t32;
		T t45 = 1.0/t33;
		T t46 = t34+1.0/2.0;
		T t47 = t35+1.0/2.0;
		T t48 = t36+1.0/2.0;
		T t49 = t37+1.0/2.0;
		T t50 = t38+1.0/2.0;
		T t51 = t46*udg1;
		T t52 = t47*udg2;
		T t53 = t48*udg3;
		T t54 = t49*udg4;
		T t55 = t50*udg5;
		T t56 = t3*t21*t41;
		T t57 = t4*t21*t42;
		T t58 = t5*t21*t43;
		T t59 = t6*t21*t44;
		T t60 = t7*t21*t45;
		T t61 = t39+t51+1.0/2.0;
		T t62 = t39+t52+1.0/2.0;
		T t63 = t39+t53+1.0/2.0;
		T t64 = t39+t54+1.0/2.0;
		T t65 = t39+t55+1.0/2.0;
		T t66 = t46+t56;
		T t67 = t47+t57;
		T t68 = t48+t58;
		T t69 = t49+t59;
		T t70 = t50+t60;
		T t81 = t40+t51+t52+t53+t54+t55+5.0/2.0;
		T t71 = t66*udg9;
		T t72 = t67*udg10;
		T t73 = t68*udg11;
		T t74 = t66*udg17;
		T t75 = t69*udg12;
		T t76 = t67*udg18;
		T t77 = t70*udg13;
		T t78 = t68*udg19;
		T t79 = t69*udg20;
		T t80 = t70*udg21;
		T t82 = 1.0/t81;
		T t83 = t82*t82;
		T t84 = t82*udg8;
		T t85 = t82*uinf1;
		T t88 = t71+t72+t73+t75+t77;
		T t89 = t74+t76+t78+t79+t80;
		T t86 = (t14*t83)/2.0;
		T t87 = (t15*t83)/2.0;
		T t90 = t84+t85;
		T t92 = t82*t88*udg6;
		T t93 = t82*t89*udg7;
		T t91 = t86+t87;
		T t94 = -t92;
		T t95 = -t93;
		T t98 = -t82*(t92-udg14);
		T t99 = -t82*(t93-udg23);
		T t101 = -odg1*(t82*(t92-udg14)+t82*(t93-udg23));
		T t96 = t94+udg14;
		T t97 = t95+udg23;
		T t100 = t98+t99;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t61*t82*udg6+nlg2*t61*t82*udg7;
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t62*t82*udg6+nlg2*t62*t82*udg7;
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t63*t82*udg6+nlg2*t63*t82*udg7;
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t64*t82*udg6+nlg2*t64*t82*udg7;
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t65*t82*udg6+nlg2*t65*t82*udg7;
		f[5*ng+i] = nlg1*(t101+uinf1+t14*t82)+tau6*(udg6-uhg6)+nlg2*t82*udg6*udg7;
		f[6*ng+i] = nlg2*(t101+uinf1+t15*t82)+tau7*(udg7-uhg7)+nlg1*t82*udg6*udg7;
		f[7*ng+i] = nlg1*(t90*udg6+t82*t101*udg6+odg1*param21*uinf20*(t71*uinf14+t72*uinf15+t73*uinf16+t75*uinf17+t77*uinf18-uinf19*(-udg16+t88*t91+t81*(t83*udg7*(udg15-t82*t88*udg7)-t83*udg6*(t92-udg14)))))+nlg2*(t90*udg7+t82*t101*udg7+odg1*param21*uinf20*(t74*uinf14+t76*uinf15+t78*uinf16+t79*uinf17+t80*uinf18-uinf19*(-udg24+t89*t91+t81*(t83*udg6*(udg22-t82*t89*udg6)-t83*udg7*(t93-udg23)))))+tau8*(udg8-uhg8);
	}
}

template <typename T> void opuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
    T fstate[16];
	T fwall[16];
	T Ustate[24];
	T Uwall[24];
	for (int i = 0; i <ng; i++) {

        T param2 = param[1];
		T param20 = param[19];
		T uinf1 = uinf[0];
		T uinf12 = uinf[11];
		T tau6 = tau[5];
		T tau7 = tau[6];

		T nlg1 = nlg[0*ng+i];
		T nlg2 = nlg[1*ng+i];
        T t2 = nlg1*udg[5*ng+i];
		T t3 = nlg2*udg[6*ng+i];
		T t4 = t2+t3;

        for (int ii = 0; ii<24; ii++)
        {
            Ustate[ii]  = udg[ii*ng+i];
            Uwall[ii]   = udg[ii*ng+i];
        }
        Uwall[5] = udg[5*ng+i]-(nlg1*t4)/2.0;
		Uwall[6] = udg[6*ng+i]-(nlg2*t4)/2.0;

        opuFluxKernel(i, fstate, xdg, Ustate, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
		opuFluxKernel(i, fwall,  xdg, Uwall,  odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
		
        for (int kk=0; kk<16; kk++)
		{
			fstate[kk] = 0.5 * (fstate[kk] + fwall[kk]);
		}

		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = 0.0;
		f[4*ng+i] = 0.0;
		f[5*ng+i] = nlg1*fstate[5] + nlg2*fstate[5+8] + tau6*(Ustate[5]-uhg[5*ng+i]);
		f[6*ng+i] = nlg1*fstate[6] + nlg2*fstate[6+8] + tau7*(Ustate[6]-uhg[6*ng+i]);
		f[7*ng+i] = 0.0;
	}
}

template <typename T> void opuFbou(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ib, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	if (ib == 1)
		opuFbou1(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 2)
		opuFbou2(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
	else if (ib == 3)
		opuFbou3(f, xdg, udg, odg, wdg, uhg, nlg, tau, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
}

template void opuFbou(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFbou(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, int, Mutation::Mixture *);
