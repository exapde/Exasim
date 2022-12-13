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
		T uinf1 = uinf[0];
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
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T uhg8 = uhg[7*ng+i];
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

		// mix->averageDiffusionCoeffs(D_i);
        // nondimensionalize_diffusionCoeffs(D_i, (double*)uinf, nspecies, ndim);
        // mix->speciesHOverRT(h_i);
        // nondimensionalize_enthalpies(h_i, (double*)uinf, nspecies, ndim);

        uinf1 = mix->P() / rhoe_scale;

		T t2 = 1.0/3.141592653589793;
		T t3 = udg1*1.0E+12;
		T t4 = udg2*1.0E+12;
		T t5 = udg3*1.0E+12;
		T t6 = udg4*1.0E+12;
		T t7 = udg5*1.0E+12;
		T t8 = atan(t3);
		T t9 = atan(t4);
		T t10 = atan(t5);
		T t11 = atan(t6);
		T t12 = atan(t7);
		T t13 = t2*t8;
		T t14 = t2*t9;
		T t15 = t2*t10;
		T t16 = t2*t11;
		T t17 = t2*t12;
		T t18 = t13+1.0/2.0;
		T t19 = t14+1.0/2.0;
		T t20 = t15+1.0/2.0;
		T t21 = t16+1.0/2.0;
		T t22 = t17+1.0/2.0;
		T t23 = t18*udg1;
		T t24 = t19*udg2;
		T t25 = t20*udg3;
		T t26 = t21*udg4;
		T t27 = t22*udg5;
		T t28 = t23+3.182454300088011E-13;
		T t29 = t24+3.182454300088011E-13;
		T t30 = t25+3.182454300088011E-13;
		T t31 = t26+3.182454300088011E-13;
		T t32 = t27+3.182454300088011E-13;
		T t33 = t23+t24+t25+t26+t27+1.591227150044006E-12;
		T t34 = 1.0/t33;
		T t35 = t34*udg8;
		T t36 = t34*uinf1;
		T t37 = t35+t36;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t28*t34*udg6+nlg2*t28*t34*udg7;
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t29*t34*udg6+nlg2*t29*t34*udg7;
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t30*t34*udg6+nlg2*t30*t34*udg7;
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t31*t34*udg6+nlg2*t31*t34*udg7;
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t32*t34*udg6+nlg2*t32*t34*udg7;
		f[5*ng+i] = nlg1*(uinf1+t34*(udg6*udg6))+tau6*(udg6-uhg6)+nlg2*t34*udg6*udg7;
		f[6*ng+i] = nlg2*(uinf1+t34*(udg7*udg7))+tau7*(udg7-uhg7)+nlg1*t34*udg6*udg7;
		f[7*ng+i] = tau8*(udg8-uhg8)+nlg1*t37*udg6+nlg2*t37*udg7;
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
		T uinf1 = uinf[0];
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
		T uhg1 = uhg[0*ng+i];
		T uhg2 = uhg[1*ng+i];
		T uhg3 = uhg[2*ng+i];
		T uhg4 = uhg[3*ng+i];
		T uhg5 = uhg[4*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
		T uhg8 = uhg[7*ng+i];
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

		// mix->averageDiffusionCoeffs(D_i);
        // nondimensionalize_diffusionCoeffs(D_i, (double*)uinf, nspecies, ndim);
        // mix->speciesHOverRT(h_i);
        // nondimensionalize_enthalpies(h_i, (double*)uinf, nspecies, ndim);

        uinf1 = mix->P() / rhoe_scale;
		T t2 = 1.0/3.141592653589793;
		T t3 = udg1*1.0E+12;
		T t4 = udg2*1.0E+12;
		T t5 = udg3*1.0E+12;
		T t6 = udg4*1.0E+12;
		T t7 = udg5*1.0E+12;
		T t8 = atan(t3);
		T t9 = atan(t4);
		T t10 = atan(t5);
		T t11 = atan(t6);
		T t12 = atan(t7);
		T t13 = t2*t8;
		T t14 = t2*t9;
		T t15 = t2*t10;
		T t16 = t2*t11;
		T t17 = t2*t12;
		T t18 = t13+1.0/2.0;
		T t19 = t14+1.0/2.0;
		T t20 = t15+1.0/2.0;
		T t21 = t16+1.0/2.0;
		T t22 = t17+1.0/2.0;
		T t23 = t18*udg1;
		T t24 = t19*udg2;
		T t25 = t20*udg3;
		T t26 = t21*udg4;
		T t27 = t22*udg5;
		T t28 = t23+3.182454300088011E-13;
		T t29 = t24+3.182454300088011E-13;
		T t30 = t25+3.182454300088011E-13;
		T t31 = t26+3.182454300088011E-13;
		T t32 = t27+3.182454300088011E-13;
		T t33 = t23+t24+t25+t26+t27+1.591227150044006E-12;
		T t34 = 1.0/t33;
		T t35 = t34*udg8;
		T t36 = t34*uinf1;
		T t37 = t35+t36;
		f[0*ng+i] = tau1*(udg1-uhg1)+nlg1*t28*t34*udg6+nlg2*t28*t34*udg7;
		f[1*ng+i] = tau2*(udg2-uhg2)+nlg1*t29*t34*udg6+nlg2*t29*t34*udg7;
		f[2*ng+i] = tau3*(udg3-uhg3)+nlg1*t30*t34*udg6+nlg2*t30*t34*udg7;
		f[3*ng+i] = tau4*(udg4-uhg4)+nlg1*t31*t34*udg6+nlg2*t31*t34*udg7;
		f[4*ng+i] = tau5*(udg5-uhg5)+nlg1*t32*t34*udg6+nlg2*t32*t34*udg7;
		f[5*ng+i] = nlg1*(uinf1+t34*(udg6*udg6))+tau6*(udg6-uhg6)+nlg2*t34*udg6*udg7;
		f[6*ng+i] = nlg2*(uinf1+t34*(udg7*udg7))+tau7*(udg7-uhg7)+nlg1*t34*udg6*udg7;
		f[7*ng+i] = tau8*(udg8-uhg8)+nlg1*t37*udg6+nlg2*t37*udg7;
	}
}

template <typename T> void opuFbou3(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uhg, T *nlg, T *tau, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
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
		T uinf1 = uinf[0];
		T tau6 = tau[5];
		T tau7 = tau[6];
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
        T udg8 = udg[7*ng+i];
		T uhg6 = uhg[5*ng+i];
		T uhg7 = uhg[6*ng+i];
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

		// mix->averageDiffusionCoeffs(D_i);
        // nondimensionalize_diffusionCoeffs(D_i, (double*)uinf, nspecies, ndim);
        // mix->speciesHOverRT(h_i);
        // nondimensionalize_enthalpies(h_i, (double*)uinf, nspecies, ndim);

        uinf1 = mix->P() / rhoe_scale;

		T t2 = 1.0/3.141592653589793;
		T t3 = udg1*1.0E+12;
		T t4 = udg2*1.0E+12;
		T t5 = udg3*1.0E+12;
		T t6 = udg4*1.0E+12;
		T t7 = udg5*1.0E+12;
		T t8 = atan(t3);
		T t9 = atan(t4);
		T t10 = atan(t5);
		T t11 = atan(t6);
		T t12 = atan(t7);
		T t13 = t2*t8;
		T t14 = t2*t9;
		T t15 = t2*t10;
		T t16 = t2*t11;
		T t17 = t2*t12;
		T t18 = t13+1.0/2.0;
		T t19 = t14+1.0/2.0;
		T t20 = t15+1.0/2.0;
		T t21 = t16+1.0/2.0;
		T t22 = t17+1.0/2.0;
		T t23 = t18*udg1;
		T t24 = t19*udg2;
		T t25 = t20*udg3;
		T t26 = t21*udg4;
		T t27 = t22*udg5;
		T t28 = t23+t24+t25+t26+t27+1.591227150044006E-12;
		T t29 = 1.0/t28;
		f[0*ng+i] = 0.0;
		f[1*ng+i] = 0.0;
		f[2*ng+i] = 0.0;
		f[3*ng+i] = 0.0;
		f[4*ng+i] = 0.0;
		f[5*ng+i] = nlg1*(uinf1+t29*(udg6*udg6))+tau6*(udg6-uhg6)+nlg2*t29*udg6*udg7;
		f[6*ng+i] = nlg2*(uinf1+t29*(udg7*udg7))+tau7*(udg7-uhg7)+nlg1*t29*udg6*udg7;
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
