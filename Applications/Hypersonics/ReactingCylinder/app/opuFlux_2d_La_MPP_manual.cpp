template <typename T> void opuFluxKernel(int i, T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
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

    T param1 = param[0];
    T param2 = param[1];
    T param19 = param[18];
    T param20 = param[19];
    T uinf1 = uinf[0];
    T uinf2 = uinf[1];
    T uinf3 = uinf[2];
    T uinf4 = uinf[3];
    T uinf5 = uinf[4];
    T uinf6 = uinf[5];
    T uinf7 = uinf[6];
    T uinf8 = uinf[7];
    T uinf9 = uinf[8];
    T uinf10 = uinf[9];
    T uinf11 = uinf[10];
    T uinf12 = uinf[11];
    T uinf13 = uinf[12];
    T uinf14 = uinf[13];
    T uinf15 = uinf[14];
    T uinf16 = uinf[15];
    T uinf17 = uinf[16];
    T uinf18 = uinf[17];
    T uinf19 = uinf[18];
    T udg1 = udg[0];
    T udg2 = udg[1];
    T udg3 = udg[2];
    T udg4 = udg[3];
    T udg5 = udg[4];
    T udg6 = udg[5];
    T udg7 = udg[6];
    T udg8 = udg[7];
    T udg9 = udg[8];
    T udg10 = udg[9];
    T udg11 = udg[10];
    T udg12 = udg[11];
    T udg13 = udg[12];
    T udg14 = udg[13];
    T udg15 = udg[14];
    T udg16 = udg[15];
    T udg17 = udg[16];
    T udg18 = udg[17];
    T udg19 = udg[18];
    T udg20 = udg[19];
    T udg21 = udg[20];
    T udg22 = udg[21];
    T udg23 = udg[22];
    T udg24 = udg[23];
    T odg1 = odg[0*ng+i];

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
    uinf2  = 0.0; //  D_1
    uinf3  = 0.0; //  D_2
    uinf4  = 0.0; //  D_3 
    uinf5  = 0.0; //  D_4
    uinf6  = 0.0; //  D_5
    uinf7  = 0.0; //  h_1
    uinf8  = 0.0; //  h_2
    uinf9  = 0.0; //  h_3
    uinf10 = 0.0; //  h_4
    uinf11 = 0.0; //  h_5 
    uinf12 = 0.0; //  mu
    uinf13 = 0.0; //  kappa

    uinf14 = 0.0;
    uinf15 = 0.0;
    uinf16 = 0.0;
    uinf17 = 0.0;
    uinf18 = 0.0;
    uinf19 = 0.0;

    T t2 = udg1*udg1;
    T t3 = udg2*udg2;
    T t4 = udg3*udg3;
    T t5 = udg4*udg4;
    T t6 = udg5*udg5;
    T t7 = 1.0/3.141592653589793;
    T t8 = udg1*1.0E+12;
    T t9 = udg2*1.0E+12;
    T t10 = udg3*1.0E+12;
    T t11 = udg4*1.0E+12;
    T t12 = udg5*1.0E+12;
    T t13 = atan(t8);
    T t14 = atan(t9);
    T t15 = atan(t10);
    T t16 = atan(t11);
    T t17 = atan(t12);
    T t33 = t2*1.0E+24;
    T t34 = t3*1.0E+24;
    T t35 = t4*1.0E+24;
    T t36 = t5*1.0E+24;
    T t37 = t6*1.0E+24;
    T t18 = t7*t13;
    T t19 = t7*t14;
    T t20 = t7*t15;
    T t21 = t7*t16;
    T t22 = t7*t17;
    T t38 = t33+1.0;
    T t39 = t34+1.0;
    T t40 = t35+1.0;
    T t41 = t36+1.0;
    T t42 = t37+1.0;
    T t23 = t18+1.0/2.0;
    T t24 = t19+1.0/2.0;
    T t25 = t20+1.0/2.0;
    T t26 = t21+1.0/2.0;
    T t27 = t22+1.0/2.0;
    T t43 = 1.0/t38;
    T t44 = 1.0/t39;
    T t45 = 1.0/t40;
    T t46 = 1.0/t41;
    T t47 = 1.0/t42;
    T t28 = t23*udg1;
    T t29 = t24*udg2;
    T t30 = t25*udg3;
    T t31 = t26*udg4;
    T t32 = t27*udg5;
    T t53 = t7*t8*t43;
    T t54 = t7*t9*t44;
    T t55 = t7*t10*t45;
    T t56 = t7*t11*t46;
    T t57 = t7*t12*t47;
    T t48 = t28+3.182454300088011E-13;
    T t49 = t29+3.182454300088011E-13;
    T t50 = t30+3.182454300088011E-13;
    T t51 = t31+3.182454300088011E-13;
    T t52 = t32+3.182454300088011E-13;
    T t58 = t23+t53;
    T t59 = t24+t54;
    T t60 = t25+t55;
    T t61 = t26+t56;
    T t62 = t27+t57;
    T t63 = t28+t29+t30+t31+t32+1.591227150044006E-12;
    T t64 = 1.0/t63;
    T t65 = t64*udg8;
    T t66 = t64*uinf1;
    T t67 = t64*udg6*udg7;
    T t68 = t65+t66;
    f[0] = odg1*t58*udg9+t48*t64*udg6;
    f[1] = odg1*t59*udg10+t49*t64*udg6;
    f[2] = odg1*t60*udg11+t50*t64*udg6;
    f[3] = odg1*t61*udg12+t51*t64*udg6;
    f[4] = odg1*t62*udg13+t52*t64*udg6;
    f[5] = uinf1+odg1*udg14+t64*(udg6*udg6);
    f[6] = t67+odg1*udg15;
    f[7] = odg1*udg16+t68*udg6;
    f[8] = odg1*t58*udg17+t48*t64*udg7;
    f[9] = odg1*t59*udg18+t49*t64*udg7;
    f[10] = odg1*t60*udg19+t50*t64*udg7;
    f[11] = odg1*t61*udg20+t51*t64*udg7;
    f[12] = odg1*t62*udg21+t52*t64*udg7;
    f[13] = t67+odg1*udg22;
    f[14] = uinf1+odg1*udg23+t64*(udg7*udg7);
    f[15] = odg1*udg24+t68*udg7;
}

template <typename T> void opuFlux(T* __restrict__ f, T* __restrict__ xdg, T* __restrict__ udg, T*__restrict__ odg, T*__restrict__ wdg, T*__restrict__ uinf, T*__restrict__ param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
    T ftmp[16];
    T utmp[16];

	for (int i = 0; i <ng; i++) {

        for (int jj = 0; jj<24; jj++)
        {
            utmp[jj] = udg[jj*ng+i];
        }

        opuFluxKernel(i, ftmp, xdg, utmp, odg, wdg, uinf, param, time, modelnumber, ng, nc, ncu, nd, ncx, nco, ncw, mix);
		
        for (int kk; kk < 16; kk++)
        {
            f[kk*ng+i] = ftmp[kk];
        }
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
