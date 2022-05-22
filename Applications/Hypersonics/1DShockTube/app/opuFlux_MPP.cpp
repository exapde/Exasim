template <typename T> void opuFlux(T *f, T *xdg, T *udg, T *odg, T *wdg, T *uinf, T *param, T time, int modelnumber, int ng, int nc, int ncu, int nd, int ncx, int nco, int ncw, Mutation::Mixture *mix)
{
	for (int i = 0; i <ng; i++) {
		T udg1 = udg[0*ng+i];
		T udg2 = udg[1*ng+i];
		T udg3 = udg[2*ng+i];
		T udg4 = udg[3*ng+i];
		T udg5 = udg[4*ng+i];
		T udg6 = udg[5*ng+i];
		T udg7 = udg[6*ng+i];
		
		int nspecies = 5;
		double rho_scale = uinf[0];
		double u_scale = uinf[1];
		double rhoe_scale = uinf[2];

		double Ucons[7] = {udg1, udg2, udg3, udg4, udg5, udg6, udg7};
		double Ustate[6];
		dimensionalizeConsVars(Ucons, (double*)uinf, nspecies, 1);
		conservativeToState(Ucons, Ustate, (double*)uinf, nspecies);
		double rhovec[5] = {Ustate[0],Ustate[1],Ustate[2],Ustate[3],Ustate[4]};
		double rhoe = Ustate[nspecies];

		mix->setState(rhovec, &rhoe, 0);
		T wdg6 = mix->P()/rhoe_scale;

		T t2 = 1.0/3.141592653589793;
		T t3 = udg1*1.0E+3;
		T t4 = udg2*1.0E+3;
		T t5 = udg3*1.0E+3;
		T t6 = udg4*1.0E+3;
		T t7 = udg5*1.0E+3;
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
		T t28 = t23+t24+t25+t26+t27+1.591548900402584E-3;
		T t29 = 1.0/t28;
		f[0*ng+i] = t29*udg6*(t23+3.183097800805168E-4);
		f[1*ng+i] = t29*udg6*(t24+3.183097800805168E-4);
		f[2*ng+i] = t29*udg6*(t25+3.183097800805168E-4);
		f[3*ng+i] = t29*udg6*(t26+3.183097800805168E-4);
		f[4*ng+i] = t29*udg6*(t27+3.183097800805168E-4);
		f[5*ng+i] = wdg6+t29*(udg6*udg6);
		f[6*ng+i] = udg6*(t29*udg7+t29*wdg6);
	}
}

template void opuFlux(double *, double *, double *, double *, double *, double *, double *, double, int, int, int, int, int, int, int, int, Mutation::Mixture *);
template void opuFlux(float *, float *, float *, float *, float *, float *, float *, float, int, int, int, int, int, int, int, int, Mutation::Mixture *);
